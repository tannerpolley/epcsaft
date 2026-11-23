#include "properties/helmholtz_model.hpp"

#include "core/trace.hpp"
#include "core/linear_error.hpp"
#include "contributions/association.hpp"
#include "contributions/born.hpp"
#include "contributions/electrolyte.hpp"
#include "contributions/helmholtz.hpp"
#include "properties/helmholtz_tape.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <condition_variable>
#include <limits>
#include <iomanip>
#include <sstream>
#include <string_view>
#include <numeric>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

namespace {
std::mutex cppad_slots_mutex;
// Runtime exit notifications may outlive an AMPL image: retain their mutexes and condition.
auto& cppad_slots = *new std::array<std::mutex, CPPAD_MAX_NUM_THREADS>;
auto& cppad_slot_exit = *new std::condition_variable;
bool cppad_parallel_enabled = false;
bool cppad_in_parallel() { return cppad_parallel_enabled; }

std::size_t cppad_thread_num() {
    thread_local const std::size_t slot = [] {
        std::lock_guard<std::mutex> lock(cppad_slots_mutex);
        for (std::size_t index = 0; index < cppad_slots.size(); ++index) {
            std::unique_lock<std::mutex> reservation(cppad_slots[index], std::try_to_lock);
            if (reservation.owns_lock()) {
                // Unlock after every TLS destructor, including earlier objects that use CppAD.
                std::notify_all_at_thread_exit(cppad_slot_exit, std::move(reservation));
                return index;
            }
        }
        throw std::runtime_error("CppAD thread slots exhausted (maximum 48 live threads)");
    }();
    return slot;
}

// Configure the library's CppAD statics before any worker can record a tape.
const bool cppad_initialized = [] {
    static std::once_flag once;
    std::call_once(once, [] {
        cppad_thread_num();
        (void)epcsaft::msa_screening_atomic();
        CppAD::thread_alloc::parallel_setup(
            CPPAD_MAX_NUM_THREADS, cppad_in_parallel, cppad_thread_num);
        CppAD::parallel_ad<double>();
        cppad_parallel_enabled = true;
    });
    return true;
}();
}  // namespace

namespace epcsaft {

struct ReferenceTaylorPreparation final {
    std::vector<std::vector<double>> path;
    std::shared_ptr<Eigen::FullPivLU<Eigen::MatrixXd>> factorization;
    PropertyActionDiagnostics diagnostics;
    ReferenceTaylorError error;
    std::vector<double> site_fractions;
    PropertyActionDiagnostics factorization_diagnostics;
    double norm = 1.0;
};

struct detail::HelmholtzEvaluation final {
    AssociationWorkspace association;
    SegmentGeometry<double> geometry;
    Eigen::FullPivLU<Eigen::MatrixXd> derivative;
    Eigen::MatrixXd sensitivity, partials;
    Eigen::VectorXd site_vector;
    std::vector<double> composition, point, weights, direction, site_weights;
    PredicateSignature signature;
    std::vector<std::size_t> coordinates;
};

HelmholtzModel::~HelmholtzModel() = default;
namespace {

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
    PropertyActionStatus* status, AssociationWorkspace* scratch
) {
    const std::size_t sites_count = sites.size();
    if (sites_count == 0) { factorization = {}; return; }  // no association sites: no mass-action system to factor
    Eigen::MatrixXd local_jacobian;
    std::vector<double> local_bonded;
    auto& jacobian = scratch ? scratch->jacobian : local_jacobian;
    auto& bonded = scratch ? scratch->products.first : local_bonded;
    bonded = bonded_densities(composition, parameters, sites, delta, std::move(bonded));
    jacobian = mass_action_jacobian(composition, parameters, sites, number_density, delta, bonded, std::move(jacobian));
    factorization.compute(jacobian);
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
        std::string message = "association Jacobian is singular";
        for (std::size_t site = 0; site < parameters.site_ids.size(); ++site)
            message += (site ? ", " : "; declared association sites: ") + parameters.site_ids[site];
        throw AssociationJacobianFailure(PropertyActionStatus::AssociationSingular, message);
    }
    if (!admitted_pivot_ratio(factorization)) {
        failure("association Jacobian is ill-conditioned");
        if (status != nullptr) {
            *status = PropertyActionStatus::AssociationIllConditioned;
            return;
        }
        throw AssociationJacobianFailure(PropertyActionStatus::AssociationIllConditioned, "association Jacobian is ill-conditioned");
    }
}

namespace {

template<class Result> Result association_result(Result result, const core::AssociationObservations& observations) {
    observations.merge(result.diagnostics.site_conditioning);
    return result;
}

// association_delta at a central state, or none for a model without sites.
std::vector<double> central_strength(const Parameters& parameters, double temperature, double density,
    const std::vector<double>& composition, const detail::SegmentGeometry<double>& geometry,
    std::vector<double> delta = {}, AssociationStrengthWorkspace<double>* scratch = nullptr) {
    return parameters.site_component.empty() ? std::vector<double>{}
        : association_delta(temperature, density, composition, parameters, geometry.sigma, geometry.diameter,
            static_cast<const std::array<double, 4>*>(nullptr), std::move(delta), scratch);
}
std::vector<double> central_strength(const Parameters& parameters, double temperature, double density,
    const std::vector<double>& composition) {
    return parameters.site_component.empty() ? std::vector<double>{} : central_strength(parameters, temperature, density,
        composition, detail::admitted_segment_geometry(temperature, parameters));
}

// One central association factorization with its rank and pivot diagnostics; a failure
// is reported through `status` and the factorization is still returned.
std::unique_ptr<Eigen::FullPivLU<Eigen::MatrixXd>> association_factorization(
    const Parameters& parameters, double density, const std::vector<double>& sites, const std::vector<double>& strength,
    const std::vector<double>& composition, PropertyActionStatus& status, PropertyActionDiagnostics& diagnostics, const policy::SolverOptions& options
) {
    if (sites.empty()) return nullptr;
    auto result = std::make_unique<Eigen::FullPivLU<Eigen::MatrixXd>>();
    properties::site_factorization(*result,
        sites, composition, parameters, density * NUMBER_DENSITY_PER_ANGSTROM3, strength, &status);
    diagnostics.site_factorizations = 1;
    diagnostics.site_rank = static_cast<std::size_t>(result->rank());
    if (!options.nested_solve) core::retain_rank(diagnostics.site_conditioning,
        core::conditioning(*result, "implicit_site_fractions", 0, "endpoint"));
    const auto pivots = result->matrixLU().diagonal().cwiseAbs();
    if (pivots.size() != 0 && pivots.minCoeff() > 0.0) diagnostics.site_pivot_ratio = pivots.maxCoeff() / pivots.minCoeff();
    return result;
}

void require(bool condition, const char* field) {
    if (!condition) {
        throw std::invalid_argument(std::string("invalid ") + field);
    }
}

void require_molar_density(double molar_density) {
    if (!std::isfinite(molar_density) || molar_density <= 0.0)
        throw std::invalid_argument("invalid molar_density " + core::optional_number(molar_density)
            + " mol/m^3; need finite rho > 0 mol/m^3");
}

// Each law slot in double: tape replay skips the domain refusal (thermal-rule volumes: the double site solve).
void require_law_domain(double T, const Parameters& parameters, std::vector<double>* scratch = nullptr) {
    std::vector<double> local;
    auto& laws = scratch ? *scratch : local;
    laws = temperature_law_values(T, parameters, std::move(laws));
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
// EqID: residual_enthalpy; residual internal energy plus the residual pressure-volume term.
Jet residual_enthalpy(const Jet& internal, const Jet& temperature, const Jet& volume, const Jet& phi_v) {
    return jet_add(internal, jet_scale(jet_multiply(jet_multiply(temperature, volume), phi_v),
        -GAS_CONSTANT_J_PER_MOL_K));
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
    const AssociationAction* association = nullptr;

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
                &coefficient_residual, &entry.outputs, nullptr, nullptr, association);
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
    if (h) result.residual[0] = residual_enthalpy(internal, temperature, volume, phi_v);
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

bool is_caloric(PropertyObservable observable) { return static_cast<std::size_t>(observable) < CALORIC_PROPERTY_COUNT; }
bool is_ideal(PropertyObservable observable) { return is_caloric(observable) && property_family(observable) == 0; }
bool is_fugacity(PropertyObservable observable) { return observable == PropertyObservable::LogFugacityCoefficient; }
bool is_volume_partial(PropertyObservable observable) { return !is_caloric(observable) && !is_fugacity(observable); }
bool is_total(PropertyObservable observable) { return is_caloric(observable) && property_family(observable) == 2; }
std::size_t property_index(PropertyObservable observable) {
    return is_volume_partial(observable) ? 3 : is_fugacity(observable) ? 0 : property_kind(observable);
}
bool is_heat_capacity(PropertyObservable observable) {
    const auto index = property_index(observable);
    return !is_fugacity(observable) && (index == 3 || index == 4);
}
void append_active_parameter_values(const HelmholtzModel& model, std::vector<double>& values) {
    for (const auto& item : model.active_parameters())
        values.push_back(detail::family_values(model.parameters(), item.family).first[item.indices.front()]);
}

// The tape's independent point (T, V, n, theta).
std::vector<double> action_point(
    const HelmholtzModel& model, double temperature, double volume,
    const std::vector<double>& amounts, bool fixed_parameters = false, std::vector<double> point = {}
) {
    point.assign({temperature, volume});
    point.insert(point.end(), amounts.begin(), amounts.end());
    if (!fixed_parameters) append_active_parameter_values(model, point);
    return point;
}

void validate_parameter_coordinate(
    const HelmholtzModel& model, ParameterFamily family,
    const std::vector<std::size_t>& indices, bool active
) {
    const auto& parameters = model.parameters();
    if (active) detail::validate_active_screening_coordinate(parameters, family);
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
            require(parameters.charge_number[indices.front()] == 0 && density_temperature_permittivity(parameters.permittivity_form),
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
    require(directions.size() <= 3, "action order");
    std::vector<bool> seen(static_cast<std::size_t>(PropertyObservable::LogFugacityCoefficient) + 1, false);
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
        require(all_finite(direction.composition), "direction composition");
        // A direction's size carries its coordinate's units, so each sum's defect is relative to its own terms.
        double balance = 0.0, charge = 0.0, size = 0.0, charge_size = 0.0;
        for (std::size_t i = 0; i < components; ++i) {
            const double term = direction.composition[i] * model.parameters().charge_number[i];
            balance += direction.composition[i];
            charge += term;
            size += std::abs(direction.composition[i]);
            charge_size += std::abs(term);
        }
        require(std::abs(balance) <= policy::COMPOSITION_SUM_TOLERANCE.value * std::max(1.0, size),
            "direction composition balance");
        require(std::abs(charge) <= policy::COMPOSITION_SUM_TOLERANCE.value * std::max(1.0, charge_size),
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
    for (const auto& [s, i, j, k] : p.cation_mixing)
        if (x[s] == 0.0 && ((component == i && x[i] == 0.0 && x[j] > 0.0)
            || (component == j && x[j] == 0.0 && x[i] > 0.0))) return true;
    return false;
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
        const auto& direction = context.directions[id];
        return direction.temperature == 0.0 && direction.molar_density == 0.0
            && std::all_of(direction.composition.begin(), direction.composition.end(),
            [](double value) { return value == 0.0; })
            && std::all_of(direction.parameter.begin(), direction.parameter.end(),
            [](double value) { return value == 0.0; });
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
void validate_cation_mixing(const Parameters& p) {
    const auto size = p.cation_mixing.size();
    require(p.cation_pair_mixing.size() == size && p.cation_triple_mixing.size() == size, "cation mixing record");
    std::set<std::pair<int,int>> pairs;
    for (const auto& [s, i, j, k] : p.cation_mixing) {
        require(pairs.emplace(i,j).second, "duplicate cation mixing pair");
        for (const int index : {s, i, j, k}) require(index >= 0 && static_cast<std::size_t>(index) < p.component_count(), "cation mixing index");
        require(i < j && p.charge_number[s] == 0 && p.charge_number[i] > 0 && p.charge_number[j] > 0
            && p.charge_number[k] < 0 && p.electrolyte_form.declared(), "cation mixing participants");
        require(std::isfinite(p.molar_mass_kg_per_mol[s]) && p.molar_mass_kg_per_mol[s] > 0.0
            && all_finite(p.cation_pair_mixing) && all_finite(p.cation_triple_mixing), "cation mixing values");
    }
}
}  // namespace

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
    require(all_finite(parameters_.k_ij), "k_ij");
    require(all_finite(parameters_.l_ij), "l_ij");
    auto l_ij = parameters_.l_ij;  // l_ij < 1 keeps the cross diameter positive; a law owns its pair's value
    for (const auto& [target, slot, law] : parameters_.temperature_law_slot)
        if (target == static_cast<int>(TemperatureLaw::Lij)) l_ij.at(slot) = l_ij.at(slot % components * components + slot / components) = 0.0;
    for (std::size_t pair = 0; pair < pairs; ++pair)
        require(l_ij[pair] < 1.0, ("l_ij for pair " + std::to_string(pair / components) + "/" + std::to_string(pair % components)).c_str());
    require(all_finite(parameters_.k_ij_reciprocal_temperature_slope),
            "k_ij_temperature_slope");
    require(all_finite(parameters_.k_ij_reference_temperature_k),
            "k_ij_reference_temperature");
    require(all_finite(parameters_.sigma_term_amplitude), "sigma_amplitude");
    require(all_finite(parameters_.sigma_term_exponent), "sigma_exponent");
    require(all_finite(parameters_.molar_mass_kg_per_mol), "molar_mass");
    require(all_finite(parameters_.relative_permittivity), "relative_permittivity");
    require(all_finite(parameters_.permittivity_term_coefficient),
            "permittivity_coefficient");
    require(all_finite(parameters_.permittivity_reference_temperature_k),
            "permittivity_reference_temperature");
    require(all_finite(parameters_.ion_suppression_coefficient),
            "ion_suppression_coefficient");
    require(all_finite(parameters_.dielectric_volume_m3_per_mol), "dielectric_volume");
    require(all_finite(parameters_.dielectric_temperature_k), "dielectric_temperature");
    require(all_finite(parameters_.dipole_scaling_debye2), "dipole_scaling");
    require(all_finite(parameters_.polarizability_scaling_angstrom3),
            "polarizability_scaling");
    require(all_finite(parameters_.correlation_integral_parameter),
            "correlation_integral");
    require(all_finite(parameters_.association_energy_over_k),
            "association_energy_over_k");
    require(all_finite(parameters_.association_volume), "association_volume");
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
        && all_finite(parameters_.temperature_law_coefficient), "temperature_law_coefficient");
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
        require(msa_screening(parameters_.electrolyte_form) || parameters_.debye_huckel_diameter_angstrom[component] > 0.0,
                "debye_huckel_diameter");
        if (born_charging(parameters_.electrolyte_form)) {
            require(parameters_.born_diameter_angstrom[component] > 0.0,
                    "born_diameter");
        }
    }
    if (born_charging(parameters_.electrolyte_form)) {
        require(parameters_.ionic_region_relative_permittivity > 0.0,
                "ionic_region_relative_permittivity");
    }
    require(std::isfinite(parameters_.ionic_region_relative_permittivity),
            "ionic_region_relative_permittivity");
    require(std::isfinite(parameters_.born_c_shell), "born_c_shell");
    require(std::isfinite(parameters_.born_c_dielectric), "born_c_dielectric");
    if (solvent_weighted_permittivity(parameters_.permittivity_form)) {
        require(parameters_.permittivity_solvent.size() == components, "permittivity solvent size");
        require(std::count(parameters_.permittivity_solvent.begin(), parameters_.permittivity_solvent.end(), 1) > 0, "permittivity solvent set");
        for (std::size_t i = 0; i < components; ++i) {
            const int member = parameters_.permittivity_solvent[i];
            require(member == 0 || member == 1, "permittivity solvent membership");
            require(!member || (parameters_.charge_number[i] == 0 && (parameters_.relative_permittivity[i] > 0.0
                || std::find(parameters_.permittivity_term_component.begin(), parameters_.permittivity_term_component.end(), i) != parameters_.permittivity_term_component.end())), "permittivity neutral solvent");
        }
    }
    if (parameters_.permittivity_form.is(PermittivityForm::SolventVolumeWeighting)) {
        const auto& volumes = parameters_.solvent_molar_volume_m3_per_mol;
        require(volumes.size() == components && all_finite(volumes) && parameters_.permittivity_solvent_normalized, "permittivity solvent volumes");
        for (std::size_t i = 0; i < components; ++i)
            require(parameters_.permittivity_solvent[i] ? volumes[i] > 0.0 : volumes[i] == 0.0, "permittivity solvent volume");
    }
    // ComponentMixing stays defined on ions alone (Bulow 2019 Eq. 7). IonSuppression and DensityTemperature
    // average neutral-solvent records; Rueben 2024 fit DipolarPerturbation ion records to water-salt mixtures.
    if (parameters_.permittivity_form.declared()
        && !ion_only_permittivity(parameters_.permittivity_form)) {
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
    association_plan_ = classify_association(parameters_, [&](std::size_t edge) {
        return seen.count({ParameterFamily::AssociationEnergyOverK, edge})
            || seen.count({ParameterFamily::AssociationVolume, edge});
    });
}

// Residual Helmholtz energy as a Taylor series in molar density at fixed T and x; Fixed = Taylor<1> adds
// the exact derivative of every coefficient along one direction in T and x.
template <typename Fixed>
static Taylor<2, Fixed> virial_energy(const Parameters& parameters, const Fixed& temperature, const std::vector<Fixed>& composition) {
    using Series = Taylor<2, Fixed>;
    const auto [sigma, diameter] = detail::admitted_segment_geometry(temperature, parameters);
    Series density;
    density.c[1] = Fixed(1.0);
    const auto delta = parameters.ideal_gas ? std::vector<Series>{} : association_delta(
        temperature, density, composition, parameters, sigma, diameter);
    for (const auto& strength : delta)
        require(std::isfinite(scalar_lead(strength)) && scalar_lead(strength) >= 0.0, "association strength");
    const auto sites = association_density_coefficients(density, composition, parameters, delta);
    const auto terms = residual_helmholtz(temperature, density, composition, parameters, sigma, diameter,
        sites, association_products(density, composition, parameters, sites, delta));
    return terms.hard_chain + terms.dispersion + terms.association + terms.dipolar;
}

// EqID: neutral_virial_coefficients
VirialCoefficientActions HelmholtzModel::virial_coefficient_actions(double temperature, const std::vector<double>& composition,
    double temperature_direction, const std::vector<double>& composition_direction) const {
    core::Scope scope("virial_coefficient_actions");
    require(std::isfinite(temperature) && temperature > 0.0, "temperature");
    validate_composition(composition);
    if (std::any_of(parameters_.charge_number.begin(), parameters_.charge_number.end(), [](int z) { return z != 0; }))
        throw std::invalid_argument("virial coefficients require a neutral model: charged declaration is nonanalytic at zero density");
    const auto energy = virial_energy(parameters_, temperature, composition);
    if (!std::isfinite(energy.c[1])) throw std::runtime_error("nonfinite second virial coefficient");
    if (!std::isfinite(2.0 * energy.c[2])) throw std::runtime_error("nonfinite third virial coefficient");
    VirialCoefficientActions result{energy.c[1], 2.0 * energy.c[2]};
    if (temperature_direction == 0.0 && composition_direction.empty()) return result;
    require(composition_direction.empty() || composition_direction.size() == composition.size(), "composition direction size");
    require(std::isfinite(temperature_direction), "temperature direction");
    double balance = 0.0;
    for (const double value : composition_direction) { require(std::isfinite(value), "composition direction"); balance += value; }
    require(std::abs(balance) <= std::numeric_limits<double>::epsilon() * static_cast<double>(composition.size()), "composition direction sum");
    using Fixed = Taylor<1>;
    std::vector<Fixed> x(composition.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        x[i].c[0] = composition[i];
        if (!composition_direction.empty()) x[i].c[1] = composition_direction[i];
    }
    Fixed t(temperature);
    t.c[1] = temperature_direction;
    const auto dual = virial_energy(parameters_, t, x);
    result.second_action = dual.c[1].c[1];
    result.third_action = 2.0 * dual.c[2].c[1];
    if (!std::isfinite(result.second_action)) throw std::runtime_error("nonfinite second virial coefficient action");
    if (!std::isfinite(result.third_action)) throw std::runtime_error("nonfinite third virial coefficient action");
    return result;
}

VirialCoefficients HelmholtzModel::virial_coefficients(double temperature, const std::vector<double>& composition) const {
    core::Scope scope("virial_coefficients");
    const auto result = virial_coefficient_actions(temperature, composition, 0.0, {});
    return {result.second, result.third};
}

void HelmholtzModel::validate_composition(
    const std::vector<double>& composition
) const {
    if (composition.size() != parameters_.component_count())
        throw std::invalid_argument("invalid composition size " + std::to_string(composition.size())
                                    + "; expected " + std::to_string(parameters_.component_count()));
    for (std::size_t i = 0; i < composition.size(); ++i) {
        if (!std::isfinite(composition[i]) || composition[i] < 0.0)
            throw std::invalid_argument("invalid composition[" + std::to_string(i) + "] fraction " + core::optional_number(composition[i])
                + "; need finite mole fraction >= 0");
    }
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
    require(all_finite(amounts), "amounts");
    require(std::all_of(amounts.begin(), amounts.end(), [](double value) {
        return value >= 0.0;
    }), "amounts");
    double total = 0.0;
    for (double amount : amounts) total += amount;
    require(std::isfinite(total) && total > 0.0, "total amount");
    if (solvent_weighted_permittivity(parameters_.permittivity_form)
        && parameters_.permittivity_solvent_normalized && detail::charged_support(amounts, parameters_)) {
        double solvent = 0.0;
        for (std::size_t i = 0; i < amounts.size(); ++i)
            if (parameters_.permittivity_solvent[i]) solvent += amounts[i];
        require(solvent > 0.0, "permittivity has no solvent denominator");
    }
}

detail::PredicateSignature HelmholtzModel::predicate_signature(
    double temperature, const std::vector<double>& composition, bool ideal, detail::PredicateSignature result
) const {
    result.clear();
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
    DerivativeOrder order, bool isothermal_hessian, bool retain_site_sensitivity, bool rank_only, const HelmholtzDerivatives* completed
) const {
    require(order != DerivativeOrder::PressureSlope, "tape pressure slope is unavailable");
    std::lock_guard<std::mutex> lock(mutex_);
    return evaluate_locked(temperature, volume, amounts, order, isothermal_hessian, retain_site_sensitivity, rank_only, completed);
}

// Declared pathway blocks are not stationary in X (pathway_weighted_association_helmholtz). Their gradient error along one
// requested column is |g1.e0| + |g0.e1| with J0 e0 = -r0 and J0 e1 = -r1 - J1 e0 (algorithms.md): after the tangent
// sweep, the Lagrangian sweep with weights (1, -mu), mu = J0^-T g0, holds g1 - J1^T mu at the sites, and g0.e1 adds -mu.r1.
// The achieved residual allowance phi = max(residual, floor) bounds each r0 component the evaluated residual cannot
// resolve: it adds |J0^-T (g1 - J1^T mu)|_1 phi.
static double coupled_site_floor(const std::vector<double>& values, const std::vector<double>& tangent,
    const std::vector<double>& second, const std::vector<double>& weights, const Eigen::FullPivLU<Eigen::MatrixXd>& lu,
    std::size_t point_size, double site_residual) {
    const auto sites = lu.rows();
    const Eigen::VectorXd error = lu.solve(-Eigen::Map<const Eigen::VectorXd>(values.data() + 1, sites));
    const Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<2>> mixed(second.data() + 2 * point_size + 1, sites);
    double floor = 0.0;
    for (Eigen::Index a = 0; a < sites; ++a)
        floor += std::abs(mixed(a) * error(a)) + std::abs(weights[1 + a] * tangent[1 + a]);
    return floor + Eigen::VectorXd(lu.transpose().solve(mixed)).lpNorm<1>() * site_residual;
}

static std::string state_coordinates(double temperature, double density) {
    std::ostringstream out;
    out << std::setprecision(17) << "; T=" << temperature << "; molar_density=" << density;
    return out.str();
}

HelmholtzDerivatives HelmholtzModel::evaluate_locked(
    double temperature,
    double volume,
    const std::vector<double>& amounts,
    DerivativeOrder order, bool isothermal_hessian, bool retain_site_sensitivity, bool rank_only, const HelmholtzDerivatives* completed
) const {
    core::AssociationObservations observations(true);
    const auto finish = [&](auto out) { observations.merge(out.rank_records); return out; };
    HelmholtzDerivatives result;
    core::CalculationWork::energy_evaluation();
    core::Scope entry("helmholtz_evaluate_locked"), scope("evaluate");
    entry.probe("requested_order", static_cast<double>(order));
    validate_point(temperature, volume, amounts);
    const std::size_t components = parameters_.component_count();
    const bool blocks = !parameters_.site_pathway.empty() || !parameters_.site_correction.empty();
    const bool phase_gradient = order == DerivativeOrder::PhaseGradient || order == DerivativeOrder::IsothermalPhaseGradient;
    const bool phase = phase_gradient || order == DerivativeOrder::PhaseHessian || order == DerivativeOrder::IsothermalPhaseHessian;
    const bool gradient = order == DerivativeOrder::Gradient || order == DerivativeOrder::FixedParameterGradient;
    const bool fixed_state = order == DerivativeOrder::FixedParameterStateHessian;
    const bool fixed_parameters = phase || fixed_state || order == DerivativeOrder::FixedParameterHessian || order == DerivativeOrder::FixedParameterGradient;
    isothermal_hessian = isothermal_hessian && !blocks && (order == DerivativeOrder::IsothermalPhaseHessian || order == DerivativeOrder::IsothermalPhaseGradient)
        && components == 2 && active_.empty() && !parameters_.ideal_gas
        && parameters_.dipole_moment_debye.empty() && parameters_.quadrupole_moment_debye_angstrom.empty()
        && (parameters_.site_component.size() == 2 || parameters_.site_component.size() == 4)
        && std::all_of(parameters_.charge_number.begin(), parameters_.charge_number.end(), [](int z) { return z == 0; })
        && std::all_of(amounts.begin(), amounts.end(), [](double n) { return n > 0.0; });
    const std::size_t result_size = 2 + components + active_.size();
    const std::size_t point_size = fixed_parameters ? 2 + components : result_size;
    double total = 0.0;
    for (double amount : amounts) total += amount;
    const double density = total / volume;
    if (!evaluation_) evaluation_ = std::make_unique<detail::HelmholtzEvaluation>();
    auto& scratch = *evaluation_;
    auto& composition = scratch.composition;
    composition.clear();
    for (double amount : amounts) composition.push_back(amount / total);
    auto& geometry = scratch.geometry;
    geometry = detail::admitted_segment_geometry(temperature, parameters_, std::move(geometry));
    require_law_domain(temperature, parameters_, &scratch.association.strength.laws);
    core::Convergence site_record;
    auto& delta = scratch.association.delta;
    delta = central_strength(parameters_, temperature, density, composition, geometry, std::move(delta), &scratch.association.strength);
    std::optional<AssociationAction> association;
    if (!parameters_.site_component.empty())
        result.site_fractions = solve_site_fractions(temperature, density,
            composition, parameters_, geometry.sigma, geometry.diameter,
            &site_record, nullptr, &delta, false, solver_options_, &result.first_numerical_failure, nullptr, &scratch.association, &association, &association_plan_);
    else (void)packing_fraction_from_geometry(density, composition, parameters_, geometry.diameter);
    const auto& sites = result.site_fractions;
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
    auto& point = scratch.point;
    point = action_point(*this, temperature, volume, amounts, fixed_parameters, std::move(point));
    auto& signature = scratch.signature;
    signature = predicate_signature(temperature, composition, false, std::move(signature));
    scope.observe([&] {
        scope.probe("predicate_signature", std::vector<double>(signature.begin(), signature.end()));
    });
    const bool reuse = completed && !scope.on() && !completed->complete_hessian
        && completed->site_fractions.size() == sites.size()
        && (sites.empty() || std::memcmp(completed->site_fractions.data(), sites.data(), sites.size() * sizeof(double)) == 0);
    if (reuse) result.value = completed->value;
    else {
        if (ensure_tape(point, sites, signature, {}, 0.0, true)) scope.probe("recording", 1.0);
        point.insert(point.end(), sites.begin(), sites.end());
        result.value = tape_->forward_zero(point, &counters_.forward).front();
        tape_->site_values = sites;
    }
    if (!std::isfinite(result.value)) throw std::runtime_error("nonfinite residual Helmholtz value" + state_coordinates(temperature, density));
    scope.probe("residual_helmholtz_extensive", result.value);
    if (order == DerivativeOrder::Value) return finish(std::move(result));

    const std::size_t sites_count = sites.size();
    const std::size_t output_count = 1 + sites_count + (tape_->ideal_available ? 5 : 0);
    auto& factorization = scratch.derivative;
    properties::site_factorization(factorization,
        sites, composition, parameters_, density * NUMBER_DENSITY_PER_ANGSTROM3, delta, nullptr, &scratch.association);
    if (reuse) {
        result.gradient = completed->gradient;
        result.hessian = completed->hessian;
        return finish(std::move(result));
    }
    // Requested chart: phase Hessians drop fixed controls, absent species and active parameters.
    auto& coordinates = scratch.coordinates;
    coordinates.clear();
    if (phase) {
        if (order == DerivativeOrder::PhaseHessian || order == DerivativeOrder::PhaseGradient) coordinates.push_back(0);
        coordinates.push_back(1);
        for (std::size_t i = 0; i < components; ++i)
            if (amounts[i] > 0.0) coordinates.push_back(2 + i);
    } else {
        coordinates.resize(point_size);
        std::iota(coordinates.begin(), coordinates.end(), std::size_t{0});
    }
    if (rank_only && !scope.on() && isothermal_hessian && coordinates.size() <= sites_count)
        coordinates.resize(1);  // Rank reads the V gradient and the independently computed isothermal Hessian.
    const std::size_t width = coordinates.size();
    auto& weights = scratch.weights;
    weights.assign(output_count, 0.0);
    weights.front() = 1.0;
    auto& site_sensitivity = scratch.sensitivity;
    site_sensitivity.resize(sites_count, width);
    result.gradient.assign(result_size, 0.0);
    if (gradient || phase_gradient || sites_count != 0) {
        const auto raw_gradient = tape_->function.Reverse(1, weights);
        ++counters_.reverse;
        bool have_site_sensitivity = false;
        // Traced calls retain the order of the gradient and Hessian observations.
        if (isothermal_hessian && !phase_gradient && !scope.on() && width <= sites_count) {
            result.hessian = isothermal(temperature, volume, amounts, &sites, &factorization, &site_sensitivity, association ? &*association : nullptr).hessian;
            have_site_sensitivity = site_sensitivity.allFinite();
        }
        if (!have_site_sensitivity) have_site_sensitivity = [&] {
            if (!association || association->strategy != AssociationPlan::ReciprocalPair || blocks || point_size != 2 + components
                || std::any_of(parameters_.charge_number.begin(), parameters_.charge_number.end(), [](int z) { return z != 0; }))
                return false;
            using D = Taylor<1>;
            const auto a = association->sites[0], b = association->sites[1];
            const std::vector<D> sigma(geometry.sigma.begin(), geometry.sigma.end());
            const std::vector<D> diameter(geometry.diameter.begin(), geometry.diameter.end());
            std::vector<D> x(components);
            site_sensitivity.setZero();  // Isolated sites remain unbonded.
            for (std::size_t k = 0; k < width; ++k) {
                const auto coordinate = coordinates[k];
                D t(temperature), v(volume), n(total);
                t.c[1] = coordinate == 0;
                v.c[1] = coordinate == 1;
                n.c[1] = coordinate >= 2;
                const D rho = n / v;
                for (std::size_t i = 0; i < components; ++i) {
                    D ni(amounts[i]);
                    ni.c[1] = coordinate == 2 + i;
                    x[i] = ni / n;
                }
                const auto sig = coordinate == 0 ? segment_sigma(t, parameters_) : sigma;
                const auto diam = coordinate == 0 ? segment_diameter(t, parameters_, sig) : diameter;
                const auto strength = association_delta(t, rho, x, parameters_, sig, diam);
                const auto coupling = [&](std::size_t i, std::size_t j) {
                    return (rho * NUMBER_DENSITY_PER_ANGSTROM3 * x[parameters_.site_component[j]]
                        * parameters_.site_multiplicity[j] * strength[i * sites_count + j]).c[1];
                };
                association->pair_tangent(sites[a], sites[b], coupling(a, b), coupling(b, a),
                    site_sensitivity(a, k), site_sensitivity(b, k));
            }
            return site_sensitivity.allFinite();
        }();
        // Site partials of the requested columns from the cheaper of forward and reverse sweeps.
        auto& partials = scratch.partials;
        partials.resize(sites_count, width);
        if (!have_site_sensitivity && width > 1 && width <= sites_count) {
            auto& direction = scratch.direction;
            direction.assign(point.size() * width, 0.0);
            for (std::size_t k = 0; k < width; ++k) direction[coordinates[k] * width + k] = 1.0;
            const auto first = tape_->function.Forward(1, width, direction);
            ++counters_.forward;
            for (std::size_t site = 0; site < sites_count; ++site)
                for (std::size_t k = 0; k < width; ++k) partials(site, k) = first[(1 + site) * width + k];
        } else if (!have_site_sensitivity && width <= sites_count) {
            auto& direction = scratch.direction;
            direction.assign(point.size(), 0.0);
            for (std::size_t k = 0; k < width; ++k) {
                direction[coordinates[k]] = 1.0;
                const auto first = tape_->function.Forward(1, direction);
                direction[coordinates[k]] = 0.0;
                ++counters_.forward;
                for (std::size_t site = 0; site < sites_count; ++site) partials(site, k) = first[1 + site];
            }
        } else if (!have_site_sensitivity && sites_count) {
            auto& site_weight = scratch.site_weights;
            site_weight.assign(output_count, 0.0);
            for (std::size_t site = 0; site < sites_count; ++site) {
                site_weight[1 + site] = 1.0;
                const auto partial = tape_->function.Reverse(1, site_weight);
                site_weight[1 + site] = 0.0;
                ++counters_.reverse;
                for (std::size_t k = 0; k < width; ++k) partials(site, k) = partial[coordinates[k]];
            }
        }
        if (!have_site_sensitivity && sites_count != 0) {
            if (association && association->strategy == AssociationPlan::BalancedBinary) {
                for (std::size_t k = 0; k < width; ++k)
                    if (!association->tangent([&](std::size_t row) { return -partials(row, k); },
                        [&](std::size_t row, double dx) { site_sensitivity(row, k) = dx; }))
                        site_sensitivity.col(k) = factorization.solve(-partials.col(k));
            } else site_sensitivity = factorization.solve(-partials);
        }
        for (std::size_t k = 0; k < width; ++k) {
            auto& derivative = result.gradient[coordinates[k]];
            derivative = raw_gradient[coordinates[k]];
            for (std::size_t site = 0; site < sites_count; ++site)
                derivative += raw_gradient[point_size + site] * site_sensitivity(site, k);
        }
        scope.probe("gradient", result.gradient);
        if (gradient) return finish(std::move(result));
        if (sites_count != 0) {
            auto& adjoint = scratch.site_vector;
            adjoint = factorization.transpose().solve(
                Eigen::Map<const Eigen::VectorXd>(raw_gradient.data() + point_size, static_cast<Eigen::Index>(sites_count)));
            for (std::size_t site = 0; site < sites_count; ++site) weights[1 + site] = -adjoint(static_cast<Eigen::Index>(site));
        }
    }
    if (sites_count != 0 && !blocks) {
        // Floor over the requested chart only: it bounds exactly the gradient entries computed here.
        auto& diagonal = scratch.site_vector;
        diagonal.resize(sites_count);
        for (std::size_t site = 0; site < sites_count; ++site)
            diagonal(site) = total * composition[parameters_.site_component[site]] * parameters_.site_multiplicity[site] / sites[site];
        result.gradient_floor = (diagonal.asDiagonal() * site_sensitivity).cwiseAbs().colwise().sum().maxCoeff() * site_residual;
    }
    if (retain_site_sensitivity) result.site_sensitivity = site_sensitivity;
    const auto complete = [this, order, isothermal_hessian, blocks, site_residual](
        HelmholtzDerivatives& out, bool defer, const auto& coordinates, const auto& point, const auto& weights,
        const auto& site_sensitivity, const Eigen::FullPivLU<Eigen::MatrixXd>* factorization, const AssociationAction* action) {
        const auto& sites = out.site_fractions;
        const std::size_t site_count = sites.size(), point_size = point.size() - site_count, size = out.gradient.size();
        core::Scope scope("helmholtz_hessian");
        ++counters_.hessian_calls;
        if (isothermal_hessian) {
            if (out.hessian.empty())
                out.hessian = isothermal(point[0], point[1], {point.begin() + 2, point.begin() + point_size}, &sites, factorization, nullptr, action).hessian;
        } else {
            out.hessian.assign(size * size, 0.0);
            auto& direction = evaluation_->direction;
            direction.assign(point.size(), 0.0);
            for (std::size_t k = 0; k < coordinates.size(); ++k) {
                const auto column = coordinates[k];
                const bool skipped = (order == DerivativeOrder::StateHessian || order == DerivativeOrder::FixedParameterStateHessian) && column >= 2 && column < 2 + parameters_.component_count();
                if (skipped && !(blocks && site_count)) continue;
                direction[column] = 1.0;
                for (std::size_t s = 0; s < site_count; ++s) direction[point_size + s] = site_sensitivity(s, k);
                const auto tangent = tape_->function.Forward(1, direction);
                ++counters_.forward;
                const auto second = tape_->function.Reverse(2, weights);
                ++counters_.reverse;
                if (blocks && site_count)
                    out.gradient_floor = std::max(out.gradient_floor, coupled_site_floor(tape_->zero_value, tangent, second, weights, *factorization, point_size, site_residual));
                if (!skipped) for (std::size_t r = 0; r < coordinates.size(); ++r) {
                    const auto row = coordinates[r];
                    if (!defer && !site_count && k == 0) out.gradient[row] = second[2 * row];
                    double value = second[2 * row + 1];
                    for (std::size_t s = 0; s < site_count; ++s) value += site_sensitivity(s, r) * second[2 * (point_size + s) + 1];
                    out.hessian[row * size + column] = value;
                }
                direction[column] = 0.0;
            }
        }
        out.volume_second_derivative = out.hessian[size + 1];
        if (order == DerivativeOrder::Hessian || order == DerivativeOrder::FixedParameterHessian) scope.probe("hessian", out.hessian);
        scope.probe("volume_second_derivative", out.volume_second_derivative);
    };
    if (phase_gradient && !blocks) {
        result.complete_hessian = [this, complete, coordinates, point, signature, weights, site_sensitivity,
            prepared = isothermal_hessian ? std::optional{std::pair{factorization, association}} : std::nullopt](HelmholtzDerivatives& out) {
            std::lock_guard lock(mutex_);
            auto& regular = evaluation_->point;
            regular.assign(point.begin(), point.end() - out.site_fractions.size());
            ensure_tape(regular, out.site_fractions, signature, {}, 0.0, true);
            tape_->forward_zero(point, &counters_.forward);
            tape_->site_values = out.site_fractions;
            complete(out, true, coordinates, point, weights, site_sensitivity, prepared ? &prepared->first : nullptr,
                prepared && prepared->second ? &*prepared->second : nullptr);
        };
    } else complete(result, false, coordinates, point, weights, site_sensitivity, &factorization, association ? &*association : nullptr);
    scope.probe("site_fractions", result.site_fractions);
    scope.probe("gradient_floor", result.gradient_floor);
    return finish(std::move(result));
}

static std::vector<double> coefficient_errors(const Parameters& p, detail::HelmholtzTape& tape,
    const std::vector<std::vector<double>>& regular, const std::vector<double>& sites,
    const std::vector<double>& coefficients, const std::vector<double>& known,
    const Eigen::FullPivLU<Eigen::MatrixXd>& lu) {
    const std::size_t count = regular.size(), ns = sites.size(), width = regular.front().size();
    if (count > 6 || !count || width != tape.regular_count) throw std::invalid_argument("coefficient error path");
    std::vector<double> errors(count, 0.0);
    if (!ns) return errors;
    // Declared pathway blocks are not stationary in X: coefficient 1 keeps its actual sensitivity to X_0.
    const bool blocks = !p.site_pathway.empty() || !p.site_correction.empty();
    const double total = std::accumulate(regular[0].begin() + 2, regular[0].begin() + 2 + p.component_count(), 0.0);
    std::vector<double> composition(regular[0].begin() + 2, regular[0].begin() + 2 + p.component_count());
    for (double& x : composition) x /= total;
    const auto sigma = segment_sigma(regular[0][0], p), diameter = segment_diameter(regular[0][0], p, sigma);
    // Differentiate the complete coefficient equations with every X_k independent. In particular,
    // d(r_k)/dX_0 includes the changing site Jacobian while holding X_k fixed.
    using AD = CppAD::AD<double>;
    std::vector<AD> independent(coefficients.begin(), coefficients.end()), output;
    CppAD::Independent(independent);
    auto lifted = tape.function.base2ad();
    for (std::size_t k = 0; k < count; ++k) {
        std::vector<AD> argument(regular[k].begin(), regular[k].end());
        argument.insert(argument.end(), independent.begin() + k * ns, independent.begin() + (k + 1) * ns);
        const auto value = lifted.Forward(k, argument);
        output.insert(output.end(), value.begin(), value.begin() + 1 + ns);
    }
    CppAD::ADFun<double> map(independent, output);
    map.Forward(0, coefficients);
    const auto jacobian = map.Jacobian(coefficients);
    // Absolute convolution terms of X_a + X_a rho_N sum_b x_b nu_b X_b Delta_ab - 1.
    using Series = Taylor<5>;
    std::vector<Series> q(width), x(p.component_count()), site(ns);
    for (std::size_t k = 0; k < count; ++k) {
        for (std::size_t j = 0; j < width; ++j) q[j].c[k] = regular[k][j];
        for (std::size_t a = 0; a < ns; ++a) site[a].c[k] = coefficients[k * ns + a];
    }
    Series amount;
    for (std::size_t i = 0; i < x.size(); ++i) amount += q[2 + i];
    for (std::size_t i = 0; i < x.size(); ++i) x[i] = q[2 + i] / amount;
    const Series density = amount / q[1];
    const auto series_parameters = detail::make_tape_parameters<Series>(p, {}, {});
    const auto series_sigma = segment_sigma(q[0], series_parameters);
    const auto series_diameter = segment_diameter(q[0], series_parameters, series_sigma);
    const auto delta = association_delta(q[0], density, x, series_parameters, series_sigma, series_diameter);
    const auto absolute = [](Series v) { for (double& c : v.c) c = std::abs(c); return v; };
    std::vector<Series> terms(ns);
    for (std::size_t a = 0; a < ns; ++a) {
        Series bonded;
        for (std::size_t b = 0; b < ns; ++b)
            bonded += absolute(x[p.site_component[b]]) * p.site_multiplicity[b] * absolute(site[b]) * absolute(delta[a * ns + b]);
        terms[a] = absolute(site[a]) + absolute(site[a]) * absolute(density) * NUMBER_DENSITY_PER_ANGSTROM3 * bonded;
        terms[a].c[0] += 1;
    }
    const auto strength = association_delta(regular[0][0], total / regular[0][1], composition, p, sigma, diameter);
    const Eigen::MatrixXd matrix = mass_action_jacobian(composition, p, sites,
        total / regular[0][1] * NUMBER_DENSITY_PER_ANGSTROM3, strength, bonded_densities(composition, p, sites, strength));
    Eigen::VectorXd previous = Eigen::VectorXd::Zero(count * ns);
    const double epsilon = std::numeric_limits<double>::epsilon();
    for (std::size_t k = 0; k < count; ++k) {
        Eigen::VectorXd right(ns);
        const Eigen::VectorXd coefficient = Eigen::Map<const Eigen::VectorXd>(coefficients.data() + k * ns, ns);
        const Eigen::VectorXd rhs = Eigen::Map<const Eigen::VectorXd>(known.data() + k * ns, ns);
        const Eigen::VectorXd defect = k == 0 ? rhs : core::compensated_defect(matrix, coefficient, rhs);
        for (std::size_t a = 0; a < ns; ++a) {
            right[a] = std::abs(defect[a]) + epsilon * std::max(terms[a].c[k], k == 0 ? double(ns) : 0.0);
            for (std::size_t j = 0; j < k * ns; ++j)
                right[a] += std::abs(jacobian[((k * (ns + 1) + 1 + a) * count * ns) + j]) * previous[j];
        }
        Eigen::VectorXd result = Eigen::VectorXd::Zero(right.size());
        for (std::ptrdiff_t i = 0; i < right.size(); ++i)
            result += lu.solve(Eigen::VectorXd::Unit(right.size(), i)).cwiseAbs() * right[i];
        previous.segment(k * ns, ns) = result;
        // The gradient's base-site contribution is the achieved-residual floor below.
        for (std::size_t j = k == 1 && !blocks ? ns : 0; j < (k + 1) * ns; ++j)
            errors[k] += std::abs(jacobian[k * (ns + 1) * count * ns + j]) * previous[j];
        if (k == 1) {
            double base_residual = epsilon * ns;
            for (std::size_t a = 0; a < ns; ++a) base_residual = std::max(base_residual, std::abs(known[a]));
            if (!blocks) {
                for (std::size_t a = 0; a < ns; ++a)
                    errors[k] += std::abs(coefficients[ns + a] * total * composition[p.site_component[a]]
                        * p.site_multiplicity[a] / sites[a]) * base_residual;
            } else {
                // The gradient floor's achieved-residual allowance: worst sign of the X_0 sensitivity over |r_0| <= phi.
                Eigen::VectorXd sensitivity(ns);
                for (std::size_t a = 0; a < ns; ++a) sensitivity[a] = jacobian[k * (ns + 1) * count * ns + a];
                const Eigen::VectorXd adjoint = lu.transpose().solve(sensitivity);
                errors[k] += adjoint.lpNorm<1>() * base_residual;
            }
        }
        if (k == 0) {
            std::vector<double> weights(count * (ns + 1), 0.0); weights[0] = 1;
            const auto hessian = map.Hessian(coefficients, weights);
            for (std::size_t a = 0; a < ns; ++a)
                for (std::size_t b = 0; b < ns; ++b)
                    errors[0] += 0.5 * previous[a] * std::abs(hessian[a * count * ns + b]) * previous[b];
        }
        if (!std::isfinite(errors[k])) throw std::runtime_error("nonfinite association coefficient error");
    }
    return errors;
}

std::vector<double> HelmholtzModel::residual_coefficient_errors(
    const std::vector<std::vector<double>>& regular, const std::vector<double>& sites) const {
    const std::size_t count = regular.size(), ns = sites.size();
    if (count > 4 || count == 0) throw std::invalid_argument("coefficient error order");
    const std::size_t width = regular.front().size();
    for (const auto& row : regular)
        if (row.size() != width || !all_finite(row)) throw std::invalid_argument("coefficient error path");
    const auto physical = 2 + parameters_.component_count();
    if (width != physical && width != physical + active_.size()) throw std::invalid_argument("coefficient error width");
    std::vector<double> errors(count, 0.0);
    if (ns == 0) return errors;
    core::CalculationWork::energy_evaluation();
    std::lock_guard<std::mutex> lock(mutex_);
    // The tape evaluates the row's active parameters, so the admission and the site system use them too.
    const auto p = detail::make_tape_parameters<double>(parameters_, active_,
        {regular[0].begin() + 2 + static_cast<std::ptrdiff_t>(parameters_.component_count()), regular[0].end()});
    double total = std::accumulate(regular[0].begin() + 2, regular[0].begin() + 2 + p.component_count(), 0.0);
    std::vector<double> composition(regular[0].begin() + 2, regular[0].begin() + 2 + p.component_count());
    for (double& x : composition) x /= total;
    const auto sigma = segment_sigma(regular[0][0], p), diameter = segment_diameter(regular[0][0], p, sigma);
    // association_delta admits the central density and geometry before the tape is recorded or replayed.
    Eigen::FullPivLU<Eigen::MatrixXd> lu;
    properties::site_factorization(lu, sites, composition, p,
        total / regular[0][1] * NUMBER_DENSITY_PER_ANGSTROM3,
        association_delta(regular[0][0], total / regular[0][1], composition, p, sigma, diameter));
    ensure_tape(regular[0], sites, predicate_signature(regular[0][0], composition), {}, 0.0, true);
    std::vector<double> coefficients(count * ns, 0.0), input(width + ns), known(count * ns, 0.0);
    std::copy(sites.begin(), sites.end(), coefficients.begin());
    for (std::size_t k = 0; k < count; ++k) {
        std::copy(regular[k].begin(), regular[k].end(), input.begin());
        std::fill(input.begin() + width, input.end(), 0.0);
        if (k == 0) std::copy(sites.begin(), sites.end(), input.begin() + width);
        auto value = k == 0 ? tape_->forward_zero(input) : tape_->function.Forward(k, input);
        for (std::size_t a = 0; a < ns; ++a) known[k * ns + a] = -value[1 + a];
        if (k == 0) continue;
        const Eigen::VectorXd next = lu.solve(Eigen::Map<const Eigen::VectorXd>(known.data() + k * ns, ns));
        for (std::size_t a = 0; a < ns; ++a) input[width + a] = coefficients[k * ns + a] = next[a];
        tape_->function.Forward(k, input);
    }
    return coefficient_errors(p, *tape_, regular, sites, coefficients, known, lu);
}

ResidualEnergyAction HelmholtzModel::residual_energy_action(double temperature, double volume,
    const std::vector<double>& amounts, const std::vector<std::vector<double>>& directions) const {
    std::lock_guard<std::mutex> lock(mutex_);
    return energy_action_locked(temperature, volume, amounts, directions);
}

ResidualEnergyGradient HelmholtzModel::repeated_direction_energy_gradient(double temperature, double volume,
    const std::vector<double>& amounts, const std::vector<double>& direction, const std::vector<std::size_t>& coordinates) const {
    std::lock_guard<std::mutex> lock(mutex_);
    ResidualEnergyGradient result;
    auto action = energy_action_locked(temperature, volume, amounts, {direction, direction}, &coordinates, &result.gradient);
    result.status = action.status;
    result.diagnostics = std::move(action.diagnostics);
    if (result.status != PropertyActionStatus::Available) result.gradient.clear();
    return result;
}

// The contraction of three to five directions, or with `coordinates` the vector D3F[v, v, e_k] of the two equal
// directions {v, v} into `gradient`.
ResidualEnergyAction HelmholtzModel::energy_action_locked(double temperature, double volume, const std::vector<double>& amounts,
    const std::vector<std::vector<double>>& directions, const std::vector<std::size_t>* coordinates, std::vector<double>* gradient) const {
    validate_point(temperature, volume, amounts);
    const std::size_t components = parameters_.component_count(), order = directions.size();
    const std::size_t width = 2 + components + active_.size();
    require(coordinates ? order == 2 : order >= 3 && order <= 5, "energy action order");
    for (const auto& direction : directions)
        require(direction.size() == width && all_finite(direction), "direction");
    double total = 0.0; for (double amount : amounts) total += amount;
    core::AssociationObservations observations(true);
    ResidualEnergyAction result;
    std::vector<ActionDirection> inserted;
    for (const auto& direction : directions)
        inserted.push_back({0.0, 0.0, {direction.begin() + 2, direction.begin() + 2 + static_cast<std::ptrdiff_t>(components)}, {}});
    if (coordinates) inserted.push_back({0.0, 0.0, std::vector<double>(components, 0.0), {}});
    for (const std::size_t k : coordinates ? *coordinates : std::vector<std::size_t>{}) {
        require(k < width, "coordinate");
        if (k >= 2 && k < 2 + components) inserted.back().composition[k - 2] = 1.0;
    }
    if ((result.status = support_status(*this, amounts, inserted)[0]) != PropertyActionStatus::Available) return association_result(std::move(result), observations);
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
    const auto refreshed = evaluate_locked(temperature, volume, amounts, DerivativeOrder::Value);
    core::retain_first(result.diagnostics.first_numerical_failure, refreshed.first_numerical_failure);
    const double density = total / volume;
    std::vector<double> composition(amounts); for (double& amount : composition) amount /= total;
    result.diagnostics.recordings = ensure_tape(point, refreshed.site_fractions,
        predicate_signature(temperature, composition), {}, 0.0, true) ? 1 : 0;
    result.status = PropertyActionStatus::Available;
    observations.drain(result.diagnostics.site_conditioning);
    const auto factorization = association_factorization(parameters_, density, refreshed.site_fractions,
        central_strength(parameters_, temperature, density, composition),
        composition, result.status, result.diagnostics, solver_options_);
    if (result.status != PropertyActionStatus::Available) return association_result(std::move(result), observations);
    if (std::any_of(norms.begin(), norms.end(), [](double norm) { return norm == 0.0; })) { result.action = 0.0;
        if (gradient) gradient->assign(coordinates->size(), 0.0); return association_result(std::move(result), observations); }
    std::optional<std::pair<std::size_t, std::size_t>> repeated;
    if (order == 3) {
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = i + 1; j < 3; ++j)
                if (normalized[i] == normalized[j] && normalized[3 - i - j] != normalized[i])
                    repeated = {{3 - i - j, i}};
    }
    double coefficient = 0.0;
    if (coordinates) {
        // EqID: repeated_direction_energy_gradient
        // D3F[v,v,.] = 2 dL/dp on the order-two path along v, L = a2 - l2.g2 - l1.g1 - l0.g0, each multiplier making L
        // stationary in the site coefficient of its order in turn (G0' l_k = dL/dX_k so far); one reverse sweep each.
        tape_->implicit_taylor({point, normalized[0], std::vector<double>(width, 0.0)}, 2, factorization.get(),
            &result.diagnostics.site_solves, &result.diagnostics.forward_sweeps, &result.diagnostics.association_residual);
        const std::size_t sites = tape_->site_count;
        std::vector<double> weight(3 * tape_->function.Range(), 0.0), adjoint(3 * (width + sites), 0.0);
        weight[2] = 1.0;
        double scale = 1.0;
        for (std::size_t sweep = 0;; ++sweep) {
            const auto partial = tape_->function.Reverse(3, weight);
            ++counters_.reverse;
            for (std::size_t i = 0; i < adjoint.size(); ++i) adjoint[i] += partial[i];
            for (std::size_t a = 0; sweep > 0 && a < sites; ++a) {  // the last multiplier's own equation
                const double residual = adjoint[3 * (width + a) + 3 - sweep];
                result.diagnostics.association_residual = std::max(result.diagnostics.association_residual,
                    std::isfinite(residual) ? std::abs(residual) / scale : std::numeric_limits<double>::infinity());
            }
            if (sweep == 3 || factorization == nullptr) break;
            Eigen::VectorXd right(sites);
            scale = 1.0;
            for (std::size_t a = 0; a < sites; ++a) scale = std::max(scale, std::abs(right[a] = adjoint[3 * (width + a) + 2 - sweep]));
            const Eigen::VectorXd multiplier = factorization->transpose().solve(right);
            ++result.diagnostics.site_solves;
            std::fill(weight.begin(), weight.end(), 0.0);
            for (std::size_t a = 0; a < sites; ++a) weight[3 * (1 + a) + 2 - sweep] = -multiplier[a];
        }
        for (const std::size_t k : *coordinates) gradient->push_back(2.0 * adjoint[3 * k] * norms[0] * norms[1]);
        ++result.diagnostics.directional_paths;
    } else if (repeated.has_value()) {
        // EqID: repeated_direction_energy_third
        // [t^5](F(p+v*t+u*t^3)-F(p+v*t-u*t^3)) = D3F[u,v,v].
        // The sparse paths avoid subtracting the singular pure-u third derivative
        // from an otherwise regular mixed response at dilute ionic states. Their difference still
        // carries D5F[v^5]/120 in each coefficient; scale v by a power of two s with s|v_i| <= n_i for
        // every present ion, keeping the path inside the radius of DH (sum n_i z_i^2)^(3/2).
        double s = 1.0;
        for (std::size_t i = 0; i < components; ++i) {
            const double v = std::abs(normalized[repeated->second][2 + i]);
            if (parameters_.charge_number[i] != 0 && amounts[i] > 0.0 && s * v > amounts[i])
                s = std::ldexp(1.0, std::ilogb(amounts[i] / v));
        }
        for (const double sign : {1.0, -1.0}) {
            std::vector<std::vector<double>> series(6, std::vector<double>(width, 0.0));
            series[0] = point;
            for (std::size_t i = 0; i < width; ++i) series[1][i] = s * normalized[repeated->second][i];
            for (std::size_t i = 0; i < width; ++i)
                series[3][i] = sign * normalized[repeated->first][i];
            coefficient += sign * tape_->implicit_taylor(series, 5, factorization.get(),
                &result.diagnostics.site_solves, &result.diagnostics.forward_sweeps,
                &result.diagnostics.association_residual).front() / (s * s);
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
    counters_.maximum_energy_order = std::max(counters_.maximum_energy_order, coordinates ? 3 : order);
    double action = repeated.has_value() ? coefficient
        : coefficient / static_cast<double>(std::size_t(1) << (order - 1));
    for (const double norm : norms) action *= norm;
    if (result.diagnostics.association_residual > policy::ASSOCIATION_ACTION_RESIDUAL.value) result.status = PropertyActionStatus::AssociationResidual;
    else if (gradient ? !all_finite(*gradient) : !std::isfinite(action)) result.status = PropertyActionStatus::Nonfinite;
    else result.action = action;
    return association_result(std::move(result), observations);
}

// Fixed-input replay uses the same contribution and stationary-association owners as the tape.
template<std::size_t D>
static Taylor<D, WiderScalar> reference_coefficients(const Parameters& p, const std::vector<ActiveParameter>& active,
    const std::vector<std::vector<double>>& path, const std::vector<double>& sites) {
    using S = Taylor<D, WiderScalar>;
    const auto ns = p.site_component.size(), width = path.front().size(), nc = p.component_count();
    std::vector<S> q(width), x(nc), xs(ns);
    for (std::size_t k = 0; k <= D; ++k) {
        for (std::size_t j = 0; j < width; ++j) q[j].c[k] = path[k][j];
        for (std::size_t j = 0; j < ns; ++j) xs[j].c[k] = sites[k * ns + j];
    }
    const auto mapped = detail::make_tape_parameters(p, active, std::vector<S>(q.begin() + 2 + nc, q.end()));
    S total;
    for (std::size_t j = 0; j < nc; ++j) total += q[2 + j];
    for (std::size_t j = 0; j < nc; ++j) x[j] = q[2 + j] / total;
    const S rho = total / q[1];
    const auto sigma = segment_sigma(q[0], mapped), diameter = segment_diameter(q[0], mapped, sigma);
    const auto delta = p.ideal_gas ? std::vector<S>{} : association_delta(q[0], rho, x, mapped, sigma, diameter);
    std::vector<double> reference(path[0].begin() + 2, path[0].begin() + 2 + nc);
    auto [electro, born] = detail::supported_electrostatics(q[0], rho, x, reference, mapped, &sigma);
    born = born_contribution(q[0], x, mapped, electro.bulk_relative_permittivity);
    return extensive_residual_helmholtz(total, residual_helmholtz(q[0], rho, x, mapped, sigma, diameter, xs,
        association_products(rho, x, mapped, xs, delta)), electro.screening, born, p);
}

static double precision_estimate(double center, WiderScalar wide, WiderScalar scale, std::size_t count) {
    if (std::numeric_limits<WiderScalar>::digits <= std::numeric_limits<double>::digits)
        throw std::runtime_error("source-reference estimate unavailable: WiderScalar is not wider than double");
    const WiderScalar nu = count * (std::numeric_limits<WiderScalar>::epsilon() / 2);
    if (!(nu < 1) || !wide_math::isfinite(wide) || !wide_math::isfinite(scale) || scale < 0)
        throw std::runtime_error("nonfinite wider reference estimate");
    const WiderScalar estimate = wide_math::abs(static_cast<WiderScalar>(center) - wide) + nu / (1 - nu) * scale
        + count * std::numeric_limits<WiderScalar>::denorm_min();
    const double error = std::nextafter(static_cast<double>(estimate), INFINITY);
    if (!wide_math::isfinite(error)) throw std::runtime_error("nonfinite wider reference estimate");
    return error;
}

template<typename Base>
static auto inverse_log_series(const std::vector<std::vector<double>>& path, std::size_t order) {
    const Base volume = path[0][1];
    std::pair<std::vector<Base>, std::vector<Base>> out;
    auto& [inverse, logarithm] = out;
    inverse.assign(order + 1, 0.0); logarithm.assign(order + 1, 0.0);
    inverse[0] = Base(1.0) / volume; logarithm[0] = wide_math::log(volume);
    for (std::size_t k = 1; k <= order; ++k) {
        logarithm[k] = Base(k) * Base(path[k][1]);
        for (std::size_t j = 1; j <= k; ++j) inverse[k] -= Base(path[j][1]) * inverse[k-j] / volume;
        for (std::size_t j = 1; j < k; ++j) logarithm[k] -= Base(j) * logarithm[j] * Base(path[k-j][1]);
        logarithm[k] /= Base(k) * volume;
    }
    return out;
}

ReferenceTaylor HelmholtzModel::reference_taylor(double temperature, double volume,
    const std::vector<double>& composition, const std::vector<double>& stoichiometry,
    const PressureActionDirection& direction, std::size_t order, bool infinite_dilution, ReferenceTaylorError* error, ReferenceTaylorBatch* batch) const {
    require(!batch || &batch->model == this, "reference batch model");
    validate_point(temperature, volume, composition);
    validate_composition(composition);
    require(order <= 4 && stoichiometry.size() == composition.size()
        && all_finite(stoichiometry) && wide_math::isfinite(direction.temperature)
        && wide_math::isfinite(direction.pressure) && direction.composition.empty(), "reference path");
    if (error && order > 1) throw std::invalid_argument("reference error order");
    core::AssociationObservations observations(true);
    ReferenceTaylor result;
    std::vector<ActionDirection> insertion{{0.0, 0.0, stoichiometry, {}}};
    require(!infinite_dilution || !detail::charged_support(composition, parameters_), "ion-free reference solvent required");
    result.status = infinite_dilution ? PropertyActionStatus::Available : support_status(*this, composition, insertion)[0];
    if (result.status != PropertyActionStatus::Available) return association_result(std::move(result), observations);
    std::lock_guard<std::mutex> lock(mutex_);
    auto key = action_point(*this, temperature, volume, composition, true);
    auto found = std::find_if(centrals_.begin(), centrals_.end(), [&](const auto& entry) { return entry.first == key; });
    if (found == centrals_.end()) {
        if (centrals_.size() == 16) centrals_.pop_front();
        core::AssociationObservations central;
        try {
            auto value = evaluate_locked(temperature, volume, composition, DerivativeOrder::FixedParameterHessian);
            central.merge(value.rank_records);
            centrals_.emplace_back(std::move(key), std::move(value));
        } catch (...) { central.forward(); throw; }
        found = std::prev(centrals_.end());
    }
    const auto& central = found->second;  // valid: centrals_ is not touched again while the lock is held
    for (const auto& rank : central.rank_records) core::AssociationObservations::publish(rank);
    const auto point = action_point(*this, temperature, volume, composition, direction.parameter.empty());
    const std::size_t width = point.size(), components = composition.size();
    std::vector<double> motion(width, 0.0);
    motion[0] = direction.temperature;
    for (const auto& parameter : direction.parameter) {
        require(wide_math::isfinite(parameter.value), "reference parameter direction");
        const auto found = std::find_if(active_.begin(), active_.end(), [&](const auto& a) {
            return a.family == parameter.family && a.indices == parameter.indices;
        });
        require(found != active_.end(), "reference active parameter");
        motion[2 + components + std::size_t(found - active_.begin())] += parameter.value;
    }
    auto path_key = point;
    path_key.insert(path_key.end(), motion.begin(), motion.end());
    path_key.insert(path_key.end(), {direction.pressure, double(order), double(infinite_dilution), double(error != nullptr)});
    std::shared_ptr<ReferenceTaylorPreparation> prepared;
    if (batch) if (const auto entry = batch->paths.find(path_key); entry != batch->paths.end())
        prepared = entry->second;
    const bool reuse = bool(prepared);
    if (!prepared) prepared = std::make_shared<ReferenceTaylorPreparation>();
    const double gas = GAS_CONSTANT_J_PER_MOL_K;
    const double pressure = gas * temperature * (1.0 / volume - central.gradient[1]);
    const double pv = gas * temperature * (-1.0 / (volume * volume) - central.hessian[central.gradient.size() + 1]);
    result.diagnostics.pressure_density_derivative = -pv * volume * volume;
    if (!wide_math::isfinite(pv) || pv >= 0.0) { result.status = PropertyActionStatus::PressureSingular; return association_result(std::move(result), observations); }
    double norm = prepared->norm;
    if (!reuse) {
        norm = wide_math::abs(direction.pressure) / (-pv * volume);
        for (std::size_t i = 0; i < width; ++i) {
            const double scale = i == 0 ? temperature : i == 1 ? volume
                : i < 2 + components ? 1.0 : std::max(1.0, wide_math::abs(point[i]));
            norm = std::max(norm, wide_math::abs(motion[i]) / scale);
        }
        if (norm == 0.0) norm = 1.0;
        prepared->norm = norm;
    }
    for (double& v : motion) v /= norm;
    const auto signature = predicate_signature(temperature, composition);
    auto& reference_tape = reference_tapes_[width != 2 + components];
    if (infinite_dilution) {
        if (!reuse) {
            if (solvent_weighted_permittivity(parameters_.permittivity_form)
                && parameters_.permittivity_solvent_normalized)
                (void)bulk_relative_permittivity(temperature, 1.0 / volume, composition, parameters_);
            const auto electro = electrolyte_contribution(temperature, 1.0 / volume, composition, composition, parameters_);
            (void)born_contribution(temperature, composition, parameters_, electro.bulk_relative_permittivity);
        }
        if (!reference_tape || reference_tape->predicates != signature || reference_tape->regular_count != width) {
            reference_tape = detail::HelmholtzTape::record(*this, point,
                central.site_fractions, signature, {}, 0.0, true, true);
            ++counters_.recordings;
            result.diagnostics.recordings = 1;
        }
    } else result.diagnostics.recordings = ensure_tape(point, central.site_fractions, signature, {}, 0.0, true) ? 1 : 0;
    auto& tape = infinite_dilution ? *reference_tape : *tape_;
    tape.site_values = central.site_fractions;
    observations.drain(result.diagnostics.site_conditioning);
    if (!reuse) {
        const auto state = batch ? batch->paths.lower_bound(point) : decltype(batch->paths)::iterator{};
        if (batch && state != batch->paths.end()
            && std::equal(point.begin(), point.end(), state->first.begin())
            && state->second->site_fractions == central.site_fractions) {
            prepared->factorization = state->second->factorization;
            const auto recordings = result.diagnostics.recordings;
            result.diagnostics = state->second->factorization_diagnostics;
            if (result.diagnostics.site_conditioning)
                result.diagnostics.site_conditioning = std::make_shared<core::NumericalRank>(*result.diagnostics.site_conditioning);
            result.diagnostics.recordings = recordings;
            result.diagnostics.site_factorizations = 0;
        } else prepared->factorization = association_factorization(parameters_, 1.0 / volume, central.site_fractions,
            central_strength(parameters_, temperature, 1.0 / volume, composition),
            composition, result.status, result.diagnostics, solver_options_);
        if (batch) {
            prepared->site_fractions = central.site_fractions;
            prepared->factorization_diagnostics = result.diagnostics;
            if (prepared->factorization_diagnostics.site_conditioning)
                prepared->factorization_diagnostics.site_conditioning = std::make_shared<core::NumericalRank>(*prepared->factorization_diagnostics.site_conditioning);
        }
    }
    const auto& factorization = prepared->factorization;
    if (result.status != PropertyActionStatus::Available) return association_result(std::move(result), observations);
    const std::size_t insertion_degree = order + 1, degree = 2 * order + 1;
    std::vector<std::vector<double>> path(degree + 1, std::vector<double>(width, 0.0));
    path[0] = point; path[1] = motion;
    // EqID: reference_single_insertion_taylor. Only one insertion reaches degree 2r+1.
    std::vector<WiderScalar> paired_wide, paired_scale;
    std::vector<double> paired_site;
    // Recurrence count for the consumed pressure (degree 3), value (1), and thermal (3) extractions.
    const std::size_t operation_count = 32 + 2 * tape.function.size_op() * 4 * (4 * 4 + 2 * 2 + 4 * 4);
    bool capture_error = false;
    const auto gradient_series = [&](const std::vector<double>& u, std::size_t last) {
        std::vector<double> gradient(last + 1, 0.0);
        if (capture_error) {
            paired_wide.assign(last + 1, 0.0L); paired_scale.assign(last + 1, 0.0L); paired_site.assign(last + 1, 0.0);
        }
        for (double sign : {1.0, -1.0}) {
            auto inserted = path;
            for (std::size_t i = 0; i < width; ++i) inserted[insertion_degree][i] += sign * u[i];
            std::vector<std::vector<double>> outputs;
            std::vector<double> sites, known;
            tape.implicit_taylor(inserted, insertion_degree + last, factorization.get(),
                &result.diagnostics.site_solves, &result.diagnostics.forward_sweeps,
                &result.diagnostics.association_residual, &outputs, capture_error ? &sites : nullptr, capture_error ? &known : nullptr);
            if (capture_error) {
                std::vector<WiderScalar> wide;
                if (degree == 1) { const auto v = reference_coefficients<1>(parameters_, active_, inserted, sites); wide.assign(v.c.begin(), v.c.end()); }
                else { const auto v = reference_coefficients<3>(parameters_, active_, inserted, sites); wide.assign(v.c.begin(), v.c.end()); }
                const auto association = central.site_fractions.empty() ? std::vector<double>(inserted.size(), 0.0)
                    : coefficient_errors(parameters_, tape, inserted, central.site_fractions, sites, known, *factorization);
                for (std::size_t k = 0; k <= last; ++k) {
                    paired_wide[k] += 0.5L * sign * wide[insertion_degree + k];
                    paired_scale[k] += 0.5L * wide_math::abs(wide[insertion_degree + k]);
                    paired_site[k] += 0.5 * association[insertion_degree + k];
                }
            }
            for (std::size_t k = 0; k <= last; ++k)
                gradient[k] += 0.5 * sign * outputs[insertion_degree + k][0];
            ++result.diagnostics.directional_paths;
        }
        return gradient;
    };
    const auto inverse_log = [&]() { return inverse_log_series<double>(path, order); };
    std::vector<double> u(width, 0.0); u[1] = volume;
    if (reuse) {
        path = prepared->path;
        const auto recordings = result.diagnostics.recordings;
        result.diagnostics = prepared->diagnostics;
        if (result.diagnostics.site_conditioning)
            result.diagnostics.site_conditioning = std::make_shared<core::NumericalRank>(*result.diagnostics.site_conditioning);
        result.diagnostics.recordings = recordings;
        result.diagnostics.site_factorizations = result.diagnostics.site_solves = 0;
        result.diagnostics.pressure_factorizations = result.diagnostics.pressure_solves = 0;
        result.diagnostics.directional_paths = result.diagnostics.forward_sweeps = 0;
        if (error) *error = prepared->error;
    } else {
        result.diagnostics.pressure_factorizations = 1;
        for (std::size_t k = 1; k <= order; ++k) {
            const auto g = gradient_series(u, k);
            const auto inverse = inverse_log().first;
            const double known = gas * (temperature * (inverse[k] - g[k] / volume)
                + motion[0] * (inverse[k-1] - g[k-1] / volume));
            path[k][1] = ((k == 1 ? direction.pressure / norm : 0.0) - known) / pv;
            ++result.diagnostics.pressure_solves;
        }
        capture_error = error != nullptr;
        const auto g = gradient_series(u, order);
        const auto [inverse, logarithm] = inverse_log();
        const auto wide_series = error ? inverse_log_series<WiderScalar>(path, order)
            : std::pair<std::vector<WiderScalar>, std::vector<WiderScalar>>{};
        if (error) {
            error->temperature_coefficient = motion[0];
            for (std::size_t k = 0; k <= order; ++k) {
                const WiderScalar r = gas, t = temperature, v = volume, t1 = motion[0];
                const auto& wi = wide_series.first;
                const WiderScalar pressure_wide = k ? r * (t * (wi[k] - paired_wide[k] / v)
                    + t1 * (wi[k-1] - paired_wide[k-1] / v)) : r * t * (wi[0] - paired_wide[0] / v);
                const WiderScalar scale = wide_math::abs(r * t) * (wide_math::abs(wi[k]) + paired_scale[k] / v)
                    + (k ? wide_math::abs(r * t1) * (wide_math::abs(wi[k-1]) + paired_scale[k-1] / v) : 0.0L);
                error->pressure[k] = k ? gas * (temperature * (inverse[k] - g[k] / volume)
                    + motion[0] * (inverse[k-1] - g[k-1] / volume)) : gas * temperature * (inverse[0] - g[0] / volume);
                error->pressure_floor[k] = precision_estimate(error->pressure[k], pressure_wide, scale, operation_count)
                    + wide_math::abs(gas * temperature / volume) * paired_site[k]
                    + (k ? wide_math::abs(gas * motion[0] / volume) * paired_site[k-1] : 0.0);
            }
        }
        for (std::size_t k = 1; k <= order; ++k) {
            const double residual = gas * (temperature * (inverse[k] - g[k] / volume)
                + motion[0] * (inverse[k-1] - g[k-1] / volume)) - (k == 1 ? direction.pressure / norm : 0.0);
            result.diagnostics.pressure_residual = std::max(result.diagnostics.pressure_residual,
                wide_math::abs(residual) / std::max(pressure, wide_math::abs(pv * path[k][1])));
        }
        observations.merge(result.diagnostics.site_conditioning);
        prepared->path = path;
        prepared->diagnostics = result.diagnostics;
        if (prepared->diagnostics.site_conditioning)
            prepared->diagnostics.site_conditioning = std::make_shared<core::NumericalRank>(*prepared->diagnostics.site_conditioning);
        if (error) prepared->error = *error;
        if (batch) batch->paths.emplace(std::move(path_key), prepared);
    }
    capture_error = error != nullptr;
    const auto [inverse, logarithm] = inverse_log();
    const auto wide_series = error ? inverse_log_series<WiderScalar>(path, order)
        : std::pair<std::vector<WiderScalar>, std::vector<WiderScalar>>{};
    u.assign(width, 0.0);
    std::copy(stoichiometry.begin(), stoichiometry.end(), u.begin() + 2);
    result.contraction = gradient_series(u, order);
    const double delta = std::accumulate(stoichiometry.begin(), stoichiometry.end(), 0.0);
    const WiderScalar wide_delta = error ? std::accumulate(stoichiometry.begin(), stoichiometry.end(), WiderScalar(0)) : WiderScalar(0);
    result.volume.resize(order + 1);
    for (std::size_t k = 0; k <= order; ++k) {
        result.contraction[k] = (result.contraction[k] - delta * logarithm[k]) * std::pow(norm, int(k));
        result.volume[k] = path[k][1] * std::pow(norm, int(k));
        if (error) {
            const WiderScalar scale = wide_math::pow(static_cast<WiderScalar>(norm), int(k));
            error->fixed[k] = precision_estimate(result.contraction[k],
                (paired_wide[k] - wide_delta * wide_series.second[k]) * scale,
                (paired_scale[k] + wide_math::abs(wide_delta * wide_series.second[k])) * scale, operation_count)
                + paired_site[k] * std::pow(norm, int(k));
        }
    }
    if (result.diagnostics.association_residual > policy::ASSOCIATION_ACTION_RESIDUAL.value)
        result.status = PropertyActionStatus::AssociationResidual;
    else if (!all_finite(result.contraction) || !all_finite(result.volume)) result.status = PropertyActionStatus::Nonfinite;
    counters_.directional_paths += result.diagnostics.directional_paths;
    counters_.forward += result.diagnostics.forward_sweeps;
    counters_.maximum_energy_order = std::max(counters_.maximum_energy_order, order + 1);
    return association_result(std::move(result), observations);
}

// ponytail: retain two signatures for alternating phases; widen only for measured multi-signature churn.
bool HelmholtzModel::ensure_tape(const std::vector<double>& point, const std::vector<double>& sites,
    const detail::PredicateSignature& signature, const std::vector<IdealCorrelation>& correlations,
    double ideal_pressure, bool residual) const {
    const auto matches = [&](const auto& tape) {
        return tape != nullptr && tape->regular_count == point.size() && tape->predicates == signature && tape->residual_available == residual
            && tape->ideal_available == !correlations.empty();
    };
    if (matches(tape_)) return false;
    if (matches(previous_tape_)) {
        tape_.swap(previous_tape_);
        return false;
    }
    previous_tape_ = detail::HelmholtzTape::record(*this, point, sites, signature, correlations, ideal_pressure, residual);
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
    core::AssociationObservations observations(true);
    PropertyActionBatch result;
    std::optional<core::Convergence> first;
    std::optional<core::Convergence> association_refusal;
    PropertyActionStatus association_status = PropertyActionStatus::Available;
    std::vector<double> sites, strength;
    std::optional<AssociationAction> association;
    const auto geometry = detail::admitted_segment_geometry(temperature, parameters_);
    (void)packing_fraction_from_geometry(molar_density, mole_fractions, parameters_, geometry.diameter);  // ideal-only actions too
    if (need_residual) {
        require_law_domain(temperature, parameters_);
        strength = central_strength(parameters_, temperature, molar_density, mole_fractions, geometry);
        try {
            sites = solve_site_fractions(temperature, molar_density, mole_fractions, parameters_, geometry.sigma,
                geometry.diameter, nullptr, nullptr, strength.empty() ? nullptr : &strength, false, solver_options_, &first, nullptr, nullptr, &association, &association_plan_);
        } catch (const core::ConvergenceFailure& error) {
            retain_density_failure(first, error, "density action association");
            association_refusal = error.record;
            need_residual = false;
            association_status = PropertyActionStatus::CentralPropertyUnavailable;
        }
    }
    const auto signature = predicate_signature(temperature, mole_fractions, ideal_available);
    core::retain_first(result.diagnostics.first_numerical_failure, first);
    result.diagnostics.recordings = ensure_tape(point, sites, signature,
        ideal_available ? correlations : std::vector<IdealCorrelation>{},
        ideal_available ? thermochemistry_->standard_pressure : 0.0, need_residual) ? 1 : 0;
    tape_->site_values = sites;
    observations.drain(result.diagnostics.site_conditioning);
    const auto factorization = need_residual ? association_factorization(parameters_, molar_density,
        tape_->site_values, strength, mole_fractions, association_status, result.diagnostics, solver_options_) : nullptr;
    ActionContext context{*tape_, point, molar_density, all,
        association_status == PropertyActionStatus::Available ? factorization.get() : nullptr};
    context.association = association ? &*association : nullptr;
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
        if (pressure_chain) {
            result.density_status = association_status;
            result.density_convergence = std::move(association_refusal);
        }
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
        if (!all_finite(result.density_action)) result.density_status = PropertyActionStatus::Nonfinite;
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
    // An empty direction list asks for values only: its order-zero jet is the value, not an action.
    for (auto& item : result.results) if (item.action) item.action = directions.empty() ? std::nullopt : std::optional<double>(*item.action * product);
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
        core::retain_first(result.diagnostics.first_numerical_failure, root->first_numerical_failure, true);
        result.diagnostics.stability_determinant = root->pressure_density_derivative
            / (GAS_CONSTANT_J_PER_MOL_K * temperature);
        result.diagnostics.pressure_residual = root->pressure_residual;
        result.diagnostics.branch = root->branch;
        result.diagnostics.returned_branch = properties::density_side(root->branch, root->molar_density, root->negative_slope_density);
        result.diagnostics.stable_root_count = root->stable_root_count;
    }
    counters_.directional_paths += context.paths.size();
    counters_.forward += context.forward_sweeps;
    counters_.maximum_energy_order = std::max(counters_.maximum_energy_order,
        result.maximum_energy_order);
    result.solver_options = solver_options_;
    return association_result(std::move(result), observations);
}

PropertyActionBatch HelmholtzModel::density_actions(
    double temperature, double molar_density, const std::vector<double>& mole_fractions,
    const std::vector<PropertyObservable>& observables,
    const std::vector<DensityActionDirection>& directions
) const {
    require(std::isfinite(temperature) && temperature > 0.0, "temperature");
    require_molar_density(molar_density);
    validate_composition(mole_fractions);
    require(!observables.empty(), "observables");
    validate_action_request(*this, observables, directions);
    std::lock_guard<std::mutex> lock(mutex_);
    return directional_actions(temperature, molar_density, mole_fractions, observables, directions);
}

// ALG-ACTIONS verification: neighbors stay on the central root through its anchor; a moved root is "branch_changed".
std::vector<core::DifferenceCheck> HelmholtzModel::pressure_action_difference(double temperature, double pressure,
    const std::vector<double>& x, DensityBranch branch, const std::vector<PropertyObservable>& observables,
    const PressureActionDirection& d, double step, double tau, double kappa) const {
    core::Scope scope("pressure_action_difference");
    core::admit_difference_parameters(step, tau, kappa);
    // Refuse unsupported directions before the exact evaluation, whose validator rejects them; shape and
    // finiteness errors still reach that validator.
    double sum = 0.0, size = 0.0;
    for (const double value : d.composition) { sum += value; size += std::abs(value); }
    std::vector<core::DifferenceCheck> out;
    if (std::any_of(d.parameter.begin(), d.parameter.end(), [](const auto& p) { return std::isfinite(p.value) && p.value != 0.0; })
        || (d.composition.size() == x.size() && std::abs(sum) > static_cast<double>(x.size()) * std::numeric_limits<double>::epsilon() * size)) {
        for (const auto o : observables) out.push_back({"observable " + std::to_string(static_cast<int>(o)), "unsupported_direction"});
        return out;
    }
    const auto exact = pressure_actions(temperature, pressure, x, branch, observables, {d});
    for (const auto& r : exact.results)
        out.push_back({"observable " + std::to_string(static_cast<int>(r.observable))
            + (r.component ? "[" + std::to_string(*r.component) + "]" : ""), "nonfinite"});
    if (!exact.molar_density) { for (auto& c : out) c.status = "neighbor_unavailable"; return out; }
    std::array<std::vector<double>, 2> values;
    for (std::size_t side = 0; side < 2; ++side) {
        const double h = side == 0 ? step : -step;
        auto w = x;
        for (std::size_t i = 0; i < d.composition.size(); ++i) w[i] += h * d.composition[i];
        const auto neighbor = pressure_actions(temperature + h * d.temperature, pressure + h * d.pressure, w, branch, observables, {d},
            exact.molar_density);
        if (!neighbor.molar_density || neighbor.results.size() != out.size()
            || std::abs(std::log(*neighbor.molar_density / *exact.molar_density)) > policy::DENSITY_ANCHOR_LOG_WINDOW.value) {
            for (auto& c : out) c.status = neighbor.molar_density ? "branch_changed" : "neighbor_unavailable";
            return out;
        }
        for (const auto& r : neighbor.results) values[side].push_back(r.value.value_or(NAN));
    }
    const double eps = std::numeric_limits<double>::epsilon();
    for (std::size_t k = 0; k < out.size(); ++k)
        out[k] = core::difference_check(out[k].label, exact.results[k].action.value_or(NAN), values[0][k], values[1][k],
            eps * std::abs(values[0][k]), eps * std::abs(values[1][k]), step, tau, kappa);
    return out;
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
    std::vector<DensityActionDirection> validation;
    validation.reserve(directions.size());
    std::vector<double> pressure_directions;
    for (const auto& direction : directions) {
        require(std::isfinite(direction.pressure), "pressure direction");
        pressure_directions.push_back(direction.pressure);
        validation.push_back({direction.temperature, 0.0, direction.composition, direction.parameter});
    }
    validate_action_request(*this, observables, validation);
    const auto unavailable = [&](const std::runtime_error& error) {
        PropertyActionBatch batch;
        batch.maximum_energy_order = directions.size() + 2;
        for (const auto observable : observables)
            for (std::size_t i = 0; i < (is_fugacity(observable) ? mole_fractions.size() : 1); ++i)
                batch.results.push_back({observable,
                    is_fugacity(observable) ? std::optional<std::size_t>(i) : std::nullopt,
                    std::nullopt, std::nullopt, PropertyActionStatus::DensityUnavailable,
                    directions.size() + (is_heat_capacity(observable) ? 2 : 1)});
        batch.density_status = PropertyActionStatus::DensityUnavailable;
        batch.solver_options = solver_options_;
        batch.density_message = error.what();
        retain_density_failure(batch.diagnostics.first_numerical_failure, error, "density closure");
        return batch;
    };
    DensityRoot root;
    try {
        root = solve_density(*this, temperature, pressure, mole_fractions, branch, anchor_density);
    } catch (const DensityStructureError& error) {
        auto batch = unavailable(error);
        batch.density_convergence = error.first_numerical_failure;
        batch.density_structure_refused = true;
        return batch;
    } catch (const core::ConvergenceFailure& error) {
        auto batch = unavailable(error);
        batch.density_convergence = error.record;
        return batch;
    } catch (const std::runtime_error& error) {
        return unavailable(error);
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
    require_molar_density(molar_density);

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
        core::CalculationWork::energy_evaluation();
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

void HelmholtzModel::reset_state_tapes() const {
    std::lock_guard<std::mutex> lock(mutex_);
    tape_.reset(); previous_tape_.reset();
}

StateBatch HelmholtzModel::states(const std::vector<double>& temperature, const std::vector<double>& input,
    const std::vector<std::vector<double>>& composition, bool at_pressure, DensityBranch branch, bool derived) const {
    core::Scope scope("states");
    const auto count = temperature.size();
    require(input.size() == 1 || input.size() == count, "state input array length");
    require(composition.size() == 1 || composition.size() == count, "composition array length");
    for (const auto& row : composition) require(row.size() == parameters_.component_count(), "composition row width");
    const bool fixed_isotherm = at_pressure && count > 1
        && has_association()
        && std::all_of(temperature.begin() + 1, temperature.end(), [&](const double& t) {
            return std::memcmp(&t, temperature.data(), sizeof(double)) == 0;
        }) && std::all_of(composition.begin() + 1, composition.end(), [&](const auto& x) {
            return x.empty() || std::memcmp(x.data(), composition.front().data(), x.size() * sizeof(double)) == 0;
        });
    std::optional<DensityScanCache> scan = fixed_isotherm ? std::optional{DensityScanCache{}} : std::nullopt;
    auto* scan_cache = scan ? &*scan : nullptr;
    StateBatch result;
    result.rows.resize(count); result.status.resize(count, StateStatus::Available); result.cause.resize(count);
    result.failure_record.resize(count); result.failure_first.resize(count);
    for (std::size_t k = 0; k < count; ++k) {
        const auto fail = [&](StateStatus status, const std::exception& error,
            std::optional<core::Convergence> record = {}, std::optional<core::Convergence> first = {}) {
            result.status[k] = status; result.cause[k] = error.what();
            result.failure_record[k] = std::move(record); result.failure_first[k] = std::move(first);
        };
        try {
            const auto& x = composition[composition.size() == 1 ? 0 : k];
            validate_electroneutrality(x);
            const double value = input[input.size() == 1 ? 0 : k];
            result.rows[k] = at_pressure ? state_at_pressure(temperature[k], value, x, branch, {}, derived,
                false, nullptr, scan_cache)
                : state(temperature[k], value, x, derived, true);
        } catch (const core::DeclarationError& e) { fail(StateStatus::DeclarationError, e, e.record);
        } catch (const core::ConvergenceFailure& e) { fail(StateStatus::ConvergenceError, e, e.record, e.first_numerical_failure);
        } catch (const DensityStructureError& e) { fail(StateStatus::DensityStructureError, e, e.first_numerical_failure, e.first_numerical_failure);
        } catch (const std::invalid_argument& e) { fail(StateStatus::ParameterError, e);
        } catch (const std::domain_error& e) { fail(StateStatus::DomainError, e);
        } catch (const std::runtime_error& e) { fail(StateStatus::EosError, e); }
    }
    return result;
}

std::vector<double> HelmholtzModel::bulk_relative_permittivity_series(const std::vector<double>& temperature,
    const std::vector<double>& density, const std::vector<std::vector<double>>& composition,
    const std::vector<double>& parameter_direction) const {
    const auto evaluate = [&](const auto& coefficient, const auto& parameters) {
        std::vector<decltype(coefficient(temperature))> x;
        for (const auto& series : composition) x.push_back(coefficient(series));
        return bulk_relative_permittivity(coefficient(temperature), coefficient(density), x, parameters);
    };
    if (temperature.size() == 1) return {evaluate([](const auto& series) { return series[0]; }, parameters_)};
    std::vector<Taylor<1>> theta;
    for (double value : active_parameter_values()) theta.emplace_back(value);
    for (std::size_t i = 0; i < parameter_direction.size(); ++i) theta[i].c[1] = parameter_direction[i];
    const auto value = evaluate([](const auto& series) {
        Taylor<1> coefficient(series[0]);
        coefficient.c[1] = series[1];
        return coefficient;
    }, detail::make_tape_parameters(parameters_, active_, theta));
    return {value.c[0], value.c[1]};
}

State HelmholtzModel::state(
    double temperature,
    double molar_density,
    const std::vector<double>& mole_fractions,
    bool derived_properties, const HelmholtzDerivatives* completed
) const {
    return state(temperature, molar_density, mole_fractions, derived_properties, false, false, completed);
}

State HelmholtzModel::state(
    double temperature,
    double molar_density,
    const std::vector<double>& mole_fractions,
    bool derived_properties, bool composition_admitted, bool fixed_parameters, const HelmholtzDerivatives* completed
) const try {
    core::AssociationObservations observations(true);
    State result;
    core::Scope scope("state");
    if (!composition_admitted) validate_composition(mole_fractions);
    require_molar_density(molar_density);
    const double volume = 1.0 / molar_density;
    const auto derivatives = [&] {
        if (parameters_.component_count() == 1 && mole_fractions[0] == 1.0 && !parameters_.ideal_gas
            && parameters_.charge_number[0] == 0.0 && parameters_.site_component.empty() && (fixed_parameters || active_.empty())) {
            validate_point(temperature, volume, mole_fractions);
            require_law_domain(temperature, parameters_);
            return pure_state_derivatives(temperature, volume);
        }
        return evaluate(temperature, volume, mole_fractions,
            scope.diagnostics_on()
                ? (fixed_parameters ? DerivativeOrder::FixedParameterHessian : DerivativeOrder::Hessian)
                : (fixed_parameters ? DerivativeOrder::FixedParameterStateHessian : DerivativeOrder::StateHessian),
            false, false, false, completed);
    }();
    const auto& sites = derivatives.site_fractions;
    const auto [electro, born, terms] = [&] {
        core::CalculationWork::energy_evaluation();
        core::Scope scope("state_residual_recompute"), evaluation("evaluate");
        const auto sigma = segment_sigma(temperature, parameters_);
        const auto diameter = segment_diameter(temperature, parameters_, sigma);
        const auto [electro, born] = detail::supported_electrostatics(
            temperature, molar_density,
            mole_fractions, mole_fractions, parameters_, &sigma
        );
        const auto terms = residual_helmholtz(
            temperature, molar_density, mole_fractions, parameters_, sigma, diameter,
            sites, association_products(molar_density, mole_fractions, parameters_, sites, parameters_.ideal_gas ? std::vector<double>{} : association_delta(temperature, molar_density, mole_fractions, parameters_, sigma, diameter)));
        return std::tuple{electro, born, terms};
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
        * (-molar_density * molar_density - phi_vv);
    const double dpressure_drho = GAS_CONSTANT_J_PER_MOL_K * temperature
        * (1.0 + (phi_vv / molar_density) / molar_density);
    const double cv_over_r = -2.0 * temperature * phi_t
        - temperature * temperature * derivatives.hessian[0];
    core::retain_first(result.first_numerical_failure, derivatives.first_numerical_failure);
    result.solver_options = solver_options_;
    result.hard_chain = terms.hard_chain;
    result.dispersion = terms.dispersion;
    result.association = terms.association;
    if (msa_screening(parameters_.electrolyte_form)) result.msa = electro.screening;
    else result.debye_huckel = electro.screening;
    result.born = born;
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
    const auto finite = [&](std::string_view field, double value, PropertyActionStatus cause = PropertyActionStatus::Nonfinite) {
        if (std::isfinite(value) && cause == PropertyActionStatus::Nonfinite) return std::optional(value);
        result.unavailable_properties.emplace_back(field, cause);
        return std::optional<double>{};
    };
    if (!result.log_fugacity_coefficient) result.unavailable_properties.emplace_back("log_fugacity_coefficient", PropertyActionStatus::CentralPropertyUnavailable);
    const auto key = [](const char* field, std::size_t index) { return std::string(field) + "[" + std::to_string(index) + "]"; };
    const std::size_t columns = fixed_parameters ? 0 : active_.size(), width = derivatives.gradient.size();
    result.density_parameter_derivative.resize(columns);
    result.log_fugacity_parameter_derivative.resize(mole_fractions.size() * columns);
    if (parameters_.permittivity_form.declared()) result.permittivity_parameter_derivative.resize(columns);
    const auto closure = !result.log_fugacity_coefficient ? PropertyActionStatus::CentralPropertyUnavailable
        : !(dpressure_drho > 0.0 && std::isfinite(dpressure_drho)) ? PropertyActionStatus::PressureSingular : PropertyActionStatus::Available;
    for (std::size_t a = 0; a < columns; ++a) {
        const std::size_t index = 2 + mole_fractions.size() + a;
        // EqID: log_fugacity_coefficient_parameter_derivative
        const double va = closure == PropertyActionStatus::Available
            ? -derivatives.hessian[width + index] / (molar_density * molar_density + phi_vv) : NAN;
        auto& density = result.density_parameter_derivative[a];
        density = finite(key("density_parameter_derivative", a), -molar_density * molar_density * va,
            closure == PropertyActionStatus::Available ? PropertyActionStatus::Nonfinite : closure);
        for (std::size_t i = 0; i < mole_fractions.size(); ++i) {
            const auto cause = result.component_property_status[i] != PropertyActionStatus::Available
                ? result.component_property_status[i] : closure != PropertyActionStatus::Available ? closure : PropertyActionStatus::Nonfinite;
            result.log_fugacity_parameter_derivative[i * columns + a] = finite(key("log_fugacity_parameter_derivative", i * columns + a),
                density && result.residual_chemical_potential_over_rt[i] ? derivatives.hessian[(2 + i) * width + index]
                    + (derivatives.hessian[(2 + i) * width + 1] - 1.0 / volume) * va : NAN, cause);
        }
    }
    if (columns && parameters_.permittivity_form.declared()) {
            std::vector<detail::AD> input{molar_density};
            const auto values = active_parameter_values();
            input.insert(input.end(), values.begin(), values.end());
            CppAD::Independent(input);
            detail::RecordingOwner recording;
            const auto record = detail::make_tape_parameters(parameters_, active_, {input.begin() + 1, input.end()});
            const auto eps = bulk_relative_permittivity(detail::AD(temperature), input[0],
                std::vector<detail::AD>(mole_fractions.begin(), mole_fractions.end()), record);
            CppAD::ADFun<double> tape(input, std::vector<detail::AD>{eps});
            recording.active = false;
            std::vector<double> point{molar_density};
            point.insert(point.end(), values.begin(), values.end());
            const auto jacobian = tape.Jacobian(point);
        for (std::size_t a = 0; a < columns; ++a) {
            const auto& density = result.density_parameter_derivative[a];
            result.permittivity_parameter_derivative[a] = finite(key("permittivity_parameter_derivative", a),
                density ? jacobian[1 + a] + jacobian[0] * *density : NAN,
                closure == PropertyActionStatus::Available ? PropertyActionStatus::Nonfinite : closure);
        }
    }
    // EqID: residual_entropy
    result.residual_entropy = GAS_CONSTANT_J_PER_MOL_K * (-derivatives.value - temperature * phi_t);
    // EqID: residual_internal_energy
    result.residual_internal_energy = -GAS_CONSTANT_J_PER_MOL_K * temperature * temperature * phi_t;
    // EqID: residual_enthalpy
    result.residual_enthalpy = residual_enthalpy(
        {result.residual_internal_energy}, {temperature}, {volume}, {phi_v}).front();
    // EqID: residual_heat_capacity_cv
    result.residual_isochoric_heat_capacity = GAS_CONSTANT_J_PER_MOL_K * cv_over_r;
    if (std::isfinite(dpressure_drho) && dpressure_drho > 0.0) {
        // EqID: residual_heat_capacity_cp
        const double cp = result.residual_isochoric_heat_capacity
            - temperature * dpressure_dt * dpressure_dt
                / dpressure_dv - GAS_CONSTANT_J_PER_MOL_K;
        result.residual_isobaric_heat_capacity = finite("residual_isobaric_heat_capacity", cp);
    } else result.unavailable_properties.emplace_back("residual_isobaric_heat_capacity", PropertyActionStatus::PressureSingular);
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
    if (!parameters_.pathway_weight.empty()) result.pathway_weight = pathway_weights(temperature, mole_fractions, parameters_)[0];
    const auto mandatory = [&](std::string_view field, double value) {
        if (!std::isfinite(value)) throw std::runtime_error("nonfinite mandatory state property: " + std::string(field) + state_coordinates(temperature, molar_density));
    };
    for (const auto& [field, value] : std::initializer_list<std::pair<const char*, double>>{
        {"molar_density", result.molar_density},
        {"hard_chain", result.hard_chain},
        {"dispersion", result.dispersion},
        {"association", result.association},
        {"debye_huckel", result.debye_huckel},
        {"msa", result.msa},
        {"born", result.born},
        {"dipolar", result.dipolar},
        {"cation_mixing", result.cation_mixing},
        {"residual_helmholtz", result.residual_helmholtz},
        {"compressibility_factor", result.compressibility_factor},
        {"pressure", result.pressure},
        {"residual_entropy", result.residual_entropy},
        {"residual_internal_energy", result.residual_internal_energy},
        {"residual_enthalpy", result.residual_enthalpy},
        {"residual_isochoric_heat_capacity", result.residual_isochoric_heat_capacity},
        {"pressure_density_derivative", result.pressure_density_derivative}}) mandatory(field, value);
    if (parameters_.permittivity_form.declared()) mandatory("bulk_relative_permittivity", result.bulk_relative_permittivity);
    for (std::size_t i = 0; i < sites.size(); ++i) mandatory(key("site_fractions", i), sites[i]);
    for (std::size_t i = 0; i < mole_fractions.size(); ++i) if (result.component_property_status[i] == PropertyActionStatus::Available)
        mandatory(key("residual_chemical_potential_over_rt", i), *result.residual_chemical_potential_over_rt[i]);
    for (const auto& [field, value] : std::initializer_list<std::pair<const char*, std::optional<double>*>>{
        {"ideal_enthalpy", &result.ideal_enthalpy},
        {"ideal_internal_energy", &result.ideal_internal_energy},
        {"ideal_entropy", &result.ideal_entropy},
        {"ideal_isobaric_heat_capacity", &result.ideal_isobaric_heat_capacity},
        {"ideal_isochoric_heat_capacity", &result.ideal_isochoric_heat_capacity},
        {"total_enthalpy", &result.total_enthalpy},
        {"total_internal_energy", &result.total_internal_energy},
        {"total_entropy", &result.total_entropy},
        {"total_isobaric_heat_capacity", &result.total_isobaric_heat_capacity},
        {"total_isochoric_heat_capacity", &result.total_isochoric_heat_capacity}}) if (*value) *value = finite(field, **value);
    if (result.log_fugacity_coefficient) for (std::size_t i = 0; i < mole_fractions.size(); ++i) {
        auto& value = (*result.log_fugacity_coefficient)[i];
        if (value) value = finite(key("log_fugacity_coefficient", i), *value);
    }
    if (!derived_properties) { observations.merge(result.rank_records); return result; }
    // Fixed-composition response and partial molar properties from the T/V Hessian columns above.
    const std::size_t components = mole_fractions.size();
    const bool slope = std::isfinite(dpressure_drho) && dpressure_drho > 0.0;
    if (!slope) for (const char* field : {"isothermal_compressibility", "isobaric_expansivity", "partial_molar_volume", "partial_molar_enthalpy", "partial_molar_entropy"})
        result.unavailable_properties.emplace_back(field, PropertyActionStatus::PressureSingular);
    const double volume_t = -dpressure_dt / dpressure_dv;  // (dv/dT)_P,x
    if (slope) {
        // EqID: isothermal_compressibility
        result.isothermal_compressibility = finite("isothermal_compressibility", 1.0 / (molar_density * dpressure_drho));
        // EqID: isobaric_expansivity
        result.isobaric_expansivity = finite("isobaric_expansivity", dpressure_dt / (molar_density * dpressure_drho));
    }
    const auto& cp = result.total_isobaric_heat_capacity;  // available only with a positive slope
    const auto& cv = result.total_isochoric_heat_capacity;
    if (cp && *cp > 0.0) {
        // EqID: joule_thomson_coefficient
        result.joule_thomson_coefficient = finite("joule_thomson_coefficient", (temperature * volume_t - volume) / *cp);
        if (cv && *cv > 0.0 && result.isothermal_compressibility) {
            // EqID: heat_capacity_ratio
            result.heat_capacity_ratio = finite("heat_capacity_ratio", *cp / *cv);
            // EqID: isentropic_compressibility
            result.isentropic_compressibility = finite("isentropic_compressibility", *result.isothermal_compressibility * *cv / *cp);
            double mass = 0.0;
            bool masses = true;
            for (std::size_t i = 0; i < components; ++i) if (mole_fractions[i] > 0.0) {
                masses = masses && parameters_.molar_mass_kg_per_mol[i] > 0.0;
                mass += mole_fractions[i] * parameters_.molar_mass_kg_per_mol[i];
            }
            // EqID: speed_of_sound
            if (masses && result.heat_capacity_ratio)
                result.speed_of_sound = finite("speed_of_sound", std::sqrt(*result.heat_capacity_ratio * dpressure_drho / mass));
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
        volumes[i] = finite(key("partial_molar_volume", i), -GAS_CONSTANT_J_PER_MOL_K * temperature * (1.0 / volume - phi_nv) / dpressure_dv);
        if (!ideal) continue;
        const auto thermal = evaluate_ideal_correlation(
            ideal_interval(thermochemistry_->component_intervals[i], temperature)->correlation,
            temperature, thermochemistry_->standard_pressure);
        // d_i = (d ln phi_i / dT)_P,x on the density chart
        const double d = phi_nt + phi_nv * volume_t + 1.0 / temperature - volume_t / volume;
        // EqID: partial_molar_enthalpy
        enthalpies[i] = finite(key("partial_molar_enthalpy", i), thermal.h - GAS_CONSTANT_J_PER_MOL_K * temperature * temperature * d);
        if (mole_fractions[i] == 0.0) continue;  // ln x_i diverges
        const double log_ideal = std::log(mole_fractions[i] * molar_density * GAS_CONSTANT_J_PER_MOL_K * temperature
            / thermochemistry_->standard_pressure) + *result.residual_chemical_potential_over_rt[i];
        // EqID: partial_molar_entropy
        entropies[i] = finite(key("partial_molar_entropy", i), thermal.s - GAS_CONSTANT_J_PER_MOL_K * (log_ideal + temperature * d));
    }
    observations.merge(result.rank_records);
    return result;
}

catch (...) { reset_state_tapes(); throw; }

// EqID: log_fugacity_coefficient_terms -- one CppAD tape of n a_k per term from the shared residual_helmholtz call.
LogFugacityTerms HelmholtzModel::log_fugacity_terms(
    double temperature, double molar_density, const std::vector<double>& mole_fractions
) const {
    core::Scope scope("log_fugacity_terms");
    using detail::AD;
    validate_electroneutrality(mole_fractions);
    require_molar_density(molar_density);
    const double volume = 1.0 / molar_density;
    validate_point(temperature, volume, mole_fractions);
    require_law_domain(temperature, parameters_);
    const auto sigma = segment_sigma(temperature, parameters_);
    const auto diameter = segment_diameter(temperature, parameters_, sigma);
    const auto sites = solve_site_fractions(temperature, molar_density, mole_fractions, parameters_, sigma, diameter,
        nullptr, nullptr, nullptr, false, solver_options_, nullptr, nullptr, nullptr, nullptr, &association_plan_);
    // Declared pathway blocks are not stationary in X: the association term also carries the implicit site chain,
    // the complete residual gradient less the fixed-site terms. Evaluated before this function records its tape.
    const bool blocks = !parameters_.site_pathway.empty() || !parameters_.site_correction.empty();
    const auto complete = blocks ? evaluate(temperature, volume, mole_fractions, DerivativeOrder::Gradient).gradient : std::vector<double>{};
    std::vector<double> point{volume};
    point.insert(point.end(), mole_fractions.begin(), mole_fractions.end());
    std::vector<AD> input(point.begin(), point.end());
    CppAD::Independent(input);
    detail::RecordingOwner recording;
    const auto parameters = detail::make_tape_parameters(parameters_, {}, std::vector<AD>{});
    const AD t(temperature);
    AD total(0.0);
    for (std::size_t i = 1; i < input.size(); ++i) total += input[i];
    const AD density = total / input[0];
    std::vector<AD> composition;
    for (std::size_t i = 1; i < input.size(); ++i) composition.push_back(input[i] / total);
    const auto sigma_ad = segment_sigma(t, parameters);
    const auto diameter_ad = segment_diameter(t, parameters, sigma_ad);
    const auto [electrolyte, born] = detail::supported_electrostatics(t, density, composition, mole_fractions, parameters, &sigma_ad);
    const std::vector<AD> site(sites.begin(), sites.end());
    const auto delta = parameters_.ideal_gas ? std::vector<AD>{} : association_delta(t, density, composition,
        parameters, sigma_ad, diameter_ad);
    const auto a = residual_helmholtz(t, density, composition, parameters, sigma_ad, diameter_ad, site,
        association_products(density, composition, parameters, site, delta));
    CppAD::ADFun<double> function(input, std::vector<AD>{total * a.hard_chain, total * a.dispersion,
        total * a.association, total * electrolyte.screening, total * born, total * a.dipolar, total * a.cation_mixing});
    recording.active = false;
    auto jacobian = function.Jacobian(point);
    const std::size_t width = point.size();
    for (std::size_t j = 0; blocks && j < width; ++j) {
        double fixed = 0.0;
        for (std::size_t k = 0; k < 7; ++k) fixed += jacobian[k * width + j];
        jacobian[2 * width + j] += complete[1 + j] - fixed;
    }
    double phi_v = 0.0;
    for (std::size_t k = 0; k < 7; ++k) phi_v += jacobian[k * width];
    const double z = 1.0 - volume * phi_v;
    require(std::isfinite(z) && z > 0.0, "compressibility factor");
    LogFugacityTerms result;
    result.log_compressibility_factor = std::log(z);
    const bool charged = detail::charged_support(mole_fractions, parameters_);
    const bool msa = msa_screening(parameters_.electrolyte_form);
    auto& screening = msa ? result.msa : result.debye_huckel;
    auto& inactive = msa ? result.debye_huckel : result.msa;
    inactive.resize(mole_fractions.size());
    for (std::size_t i = 0; i < mole_fractions.size(); ++i)
        if (!ionic_insertion_value(parameters_, mole_fractions, i, charged)) inactive[i] = 0.0;
    std::vector<std::optional<double>>* const fields[] = {&result.hard_chain, &result.dispersion,
        &result.association, &screening, &result.born, &result.dipolar, &result.cation_mixing};
    for (std::size_t k = 0; k < 7; ++k) {
        fields[k]->resize(mole_fractions.size());
        for (std::size_t i = 0; i < mole_fractions.size(); ++i)
            if (!ionic_insertion_value(parameters_, mole_fractions, i, charged)) (*fields[k])[i] = jacobian[k * width + 1 + i];
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
    bool derived_properties, bool fixed_parameters, const State* feed_state, DensityScanCache* scan
) const {
    core::Scope scope("state_at_pressure");
    const auto root = solve_density(
        *this, temperature, pressure, mole_fractions, branch, anchor_density, nullptr, scan
    );
    const bool reuse = feed_state && !scope.on() && !core::CalculationWork::current()
        && !core::AssociationObservations::current() && !derived_properties && !fixed_parameters
        && std::memcmp(&root.molar_density, &feed_state->molar_density, sizeof(double)) == 0;
    State result = reuse ? *feed_state : state(temperature, root.molar_density, mole_fractions, derived_properties, true, fixed_parameters);
    if (result.log_fugacity_coefficient) {
        // Fixed-pressure fugacity uses the prescribed P; EOS pressure retains the root residual.
        const double log_z = std::log(pressure / (GAS_CONSTANT_J_PER_MOL_K * temperature * root.molar_density));
        for (std::size_t i = 0; i < mole_fractions.size(); ++i) if (result.residual_chemical_potential_over_rt[i])
            (*result.log_fugacity_coefficient)[i] = *result.residual_chemical_potential_over_rt[i] - log_z;
    }
    result.density_convergence = root.convergence;
    for (const auto& rank : root.rank_records) core::retain_rank(result.rank_records, rank);
    core::retain_first(result.first_numerical_failure, root.first_numerical_failure, true);
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
    return admitted_packing_coefficient(temperature, mole_fractions);
}

double HelmholtzModel::admitted_packing_coefficient(double temperature, const std::vector<double>& mole_fractions) const {
    if (parameters_.ideal_gas) return 0.0;
    return packing_fraction_from_geometry(1.0, mole_fractions, parameters_,
        detail::admitted_segment_geometry(temperature, parameters_).diameter);
}

double HelmholtzModel::molar_density_upper_bound(
    double temperature,
    const std::vector<double>& mole_fractions
) const {
    core::Scope scope("molar_density_upper_bound");
    validate_composition(mole_fractions);
    require(std::isfinite(temperature) && temperature > 0.0, "temperature");
    const auto geometry = detail::admitted_segment_geometry(temperature, parameters_);
    const double coefficient = parameters_.ideal_gas ? 0.0
        : packing_fraction_from_geometry(1.0, mole_fractions, parameters_, geometry.diameter);
    return detail::density_upper_bound(coefficient, mole_fractions, geometry.sigma, parameters_);
}

std::vector<double> HelmholtzModel::active_parameter_values() const {
    std::vector<double> result;
    append_active_parameter_values(*this, result);
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
