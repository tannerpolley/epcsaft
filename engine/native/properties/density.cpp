#include "properties/helmholtz_model.hpp"
#include "pcsaft_superancillary_table.hpp"

#include "core/convergence.hpp"
#include "core/numerical_policy.hpp"
#include "core/physical_constants.hpp"
#include "core/trace.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <cstddef>
#include <iterator>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

namespace epcsaft {
namespace {

// The diagnostic scan's evaluation or time budget ran out; it ends the scan instead of
// becoming one sample's cause.
struct DiagnosticBudgetExhausted final : std::runtime_error {
    DiagnosticBudgetExhausted() : std::runtime_error("density diagnostic budget exhausted") {}
};

// 2 - phi: the golden-section interior fraction, a mathematical constant.
constexpr double GOLDEN_SECTION = 0.3819660112501051;

struct PressurePoint final {
    double residual;
    double pressure;
    double slope;
    double gradient_floor;
};

// Request-local state exists only after a resolved sampled contradiction.
struct ScanPoint final { double density; PressurePoint point; bool linked; };
struct DensityPartition final {
    static constexpr auto cap = policy::DENSITY_BRACKET_INTERVALS.value * policy::DENSITY_MAX_ITERATIONS.value;
    double lower, upper;
    std::string cause;
    std::map<double, PressurePoint> samples;
    std::size_t evaluations = 0, intervals = 0, stable = 0;
    bool active = false;
    DensityStructureError refused(const std::string& cause_now, const std::optional<core::Convergence>& first) const {
        return DensityStructureError(cause + "; enumeration interval=[" + core::optional_number(lower) + "," + core::optional_number(upper)
            + "] mol/m3; cells=" + std::to_string(intervals) + "; stable=" + std::to_string(stable)
            + "; evaluations=" + std::to_string(evaluations) + "; unresolved: " + cause_now, first);
    }

};

// P(rho) at fixed T and amounts: the shared residual Helmholtz templates as a Taylor series in volume
// at the tape's PressureSlope chart (unit total amount), P = RT (rho - A_V), dP/drho = RT (1 + A_VV/rho^2)
// and d2P/drho2 = -RT (A_VVV/rho^4 + 2 A_VV/rho^3). Only density-dependent arithmetic carries the series;
// temperature and composition terms stay double. The solved site fractions enter with their implicit volume
// coefficients from the factorized mass action, so the series is the stationary association potential's.
// It reads no tape and takes no lock.
template <typename Base = double>
class PressureSeries final {
public:
    PressureSeries(const HelmholtzModel& model, Base temperature, const std::vector<Base>& amounts)
        : model_(model), parameters_([&]() -> decltype(auto) {
            if constexpr (std::is_same_v<Base, double>) return model.parameters();
            else {
                ParameterRecord<Base> mapped;
                static_cast<FixedParameters&>(mapped) = model.parameters();
                for (std::size_t f = 0; f < detail::PARAMETER_FAMILY_COUNT; ++f) {
                    const auto family = static_cast<detail::ParameterFamily>(f);
                    const auto source = detail::family_values(model.parameters(), family);
                    const auto target = detail::family_values(mapped, family, source.second);
                    std::copy_n(source.first, source.second, target.first);
                }
                return mapped;
            }
        }()), temperature_(temperature), amounts_(std::is_same_v<Base, double> ? std::vector<Base>{} : amounts) {
        for (Base amount : amounts) total_ += amount;
        for (Base amount : amounts) composition_.push_back(amount / total_);
        sigma_ = segment_sigma(temperature, parameters_);
        if (!std::all_of(sigma_.begin(), sigma_.end(), [](Base v) { return wide_math::isfinite(v) && v > 0.0; }))
            throw std::invalid_argument("invalid sigma correlation");
        diameter_ = segment_diameter(temperature, parameters_, sigma_);
        if (!std::all_of(diameter_.begin(), diameter_.end(), [](Base v) { return wide_math::isfinite(v) && v > 0.0; }))
            throw std::invalid_argument("invalid segment diameter");
        if constexpr (!std::is_same_v<Base, double>) {
            const auto laws = temperature_law_values(temperature_, parameters_);
            const auto checked = narrow_values(laws.begin(), laws.end());
            for (const auto& [target, slot, law] : parameters_.temperature_law_slot)
                (void)temperature_law_value(static_cast<double>(temperature_), model.parameters(), checked,
                    static_cast<TemperatureLaw>(target), static_cast<std::size_t>(slot),
                    target <= int(TemperatureLaw::Lij) ? parameters_.component_count() : parameters_.site_component.size());
        }
        // With one partner per site the mass action's pair step is its exact solution at every density.
        const std::size_t sites = parameters_.site_component.size();
        for (std::size_t a = 0; a < sites; ++a)
            pairs_ = pairs_ && std::count_if(parameters_.association_active.begin() + static_cast<std::ptrdiff_t>(a * sites),
                parameters_.association_active.begin() + static_cast<std::ptrdiff_t>((a + 1) * sites), [](int v) { return v != 0; }) == 1;
    }
    mutable std::optional<DensityPartition> partition;
    mutable std::optional<core::Convergence> first_numerical_failure;
    const policy::SolverOptions& solver_options() const { return model_.solver_options(); }
    const HelmholtzModel& model() const { return model_; }

    // Coefficients of A_res/(RT) in t at V = 1/density + t, and the association gradient floor. `start` begins the
    // site solve from a neighbouring solution; `sites` receives the site fractions' coefficients in t (N per site).
    template <std::size_t N>
    Taylor<N, Base> residual(Base density, Base& floor, const std::vector<double>* start = nullptr,
                       std::vector<Base>* sites_out = nullptr, const std::vector<Base>& direction = {}) const {
        using Series = Taylor<N, Base>;
        core::Scope entry("pressure_series_residual"), scope("evaluate");
        entry.probe("requested_order", static_cast<double>(N - 1));
        const Base volume = 1.0 / density, rho0 = total_ / volume;
        Series v(volume);
        v.c[1] = direction.empty() ? Base(1) : direction[0];
        const auto path = [&] {
            if constexpr (std::is_same_v<Base, double>) return std::pair<Base, const std::vector<Base>&>(total_, composition_);
            else {
                Series total(0); std::vector<Series> composition;
                for (std::size_t i = 0; i < amounts_.size(); ++i) {
                    composition.emplace_back(amounts_[i]);
                    if (!direction.empty()) composition.back().c[1] = direction[i + 1];
                    total += composition.back();
                }
                for (auto& x : composition) x = x / total;
                return std::make_pair(total, composition);
            }
        }();
        const auto& [total, composition] = path;
        const Series rho = total / v;
        if constexpr (!std::is_same_v<Base, double>)
            if (packing_moments(rho, composition, parameters_, std::vector<Series>(diameter_.begin(), diameter_.end()))[3].c[0] >= policy::CLOSE_PACKING_FRACTION.value)
                throw std::invalid_argument("packing fraction is at or above close packing");
        // Per-thread scratch: the closure evaluates in a tight loop, so its small buffers are reused.
        static thread_local std::vector<Series> delta, site_series, residuals;
        static thread_local std::vector<double> strength;
        static thread_local Eigen::VectorXd rows, coefficient, cross, unit;
        delta.clear(); strength.clear(); site_series.clear(); residuals.clear();
        if (!parameters_.site_component.empty()) {
            if constexpr (std::is_same_v<Base, double>)
                delta = association_delta_values(temperature_, rho, composition_, parameters_, sigma_, diameter_,
                    parameters_.association_energy_over_k, parameters_.association_volume);
            else delta = association_delta_values(Series(temperature_), rho, composition, parameters_,
                std::vector<Series>(sigma_.begin(), sigma_.end()), std::vector<Series>(diameter_.begin(), diameter_.end()),
                std::vector<Series>(parameters_.association_energy_over_k.begin(), parameters_.association_energy_over_k.end()),
                std::vector<Series>(parameters_.association_volume.begin(), parameters_.association_volume.end()));
            for (const auto& d : delta) {
                if constexpr (!std::is_same_v<Base, double>)
                    if (!(d.c[0] >= 0) || !wide_math::isfinite(d.c[0])) throw std::invalid_argument("invalid extended association strength");
                strength.push_back(static_cast<double>(d.c[0]));
            }
        }
        core::Convergence record;
        const auto sites = [&]() -> std::vector<double> { try {
            if (strength.empty()) return std::vector<double>{};
            if constexpr (!std::is_same_v<Base, double>) if (!central_sites_.empty()) return narrow_values(central_sites_.begin(), central_sites_.end());
            const auto* seed = !pairs_ && start != nullptr && start->size() == parameters_.site_component.size() ? start : nullptr;
            if constexpr (std::is_same_v<Base, double>)
                return solve_site_fractions(temperature_, rho0, composition_, parameters_, sigma_, diameter_, &record, seed, &strength, true, solver_options(), &first_numerical_failure);
            else return solve_site_fractions(static_cast<double>(temperature_), static_cast<double>(rho0),
                narrow_values(composition_.begin(), composition_.end()), model_.parameters(), narrow_values(sigma_.begin(), sigma_.end()),
                narrow_values(diameter_.begin(), diameter_.end()), &record, seed, &strength, true, solver_options(), &first_numerical_failure);
        } catch (const std::exception& error) { core::retain_failure(first_numerical_failure, error, "density association"); throw; }
        }();
        site_series.clear();
        for (double site : sites) site_series.emplace_back(Base(site));
        for (std::size_t a = 0; a < central_sites_.size(); ++a) site_series[a].c[0] = central_sites_[a];
        floor = 0.0;
        if (!sites.empty()) {
            static thread_local Eigen::FullPivLU<Eigen::MatrixXd> factorization;
            if constexpr (std::is_same_v<Base, double>)
                properties::site_factorization(factorization, sites, composition_, parameters_, rho0 * NUMBER_DENSITY_PER_ANGSTROM3, strength);
            else if (central_sites_.empty()) properties::site_factorization(factorization, sites, narrow_values(composition_.begin(), composition_.end()), model_.parameters(),
                static_cast<double>(rho0 * NUMBER_DENSITY_PER_ANGSTROM3), strength);
            const auto count = static_cast<Eigen::Index>(sites.size());
            rows.resize(count); cross.resize(count);
            if constexpr (!std::is_same_v<Base, double>) {
                if (central_sites_.empty()) { central_sites_.assign(sites.begin(), sites.end()); used_ = record.iterations; }
                auto& used = used_;
                for (std::size_t k = 0; k <= N; ++k) {
                    Base normalization = 1.0;
                    for (std::size_t j = 0; direction.empty() && j < k; ++j) normalization *= volume;
                    for (;;) {
                        residuals = association_residuals(rho, composition, parameters_, site_series, delta);
                        Base norm = 0.0, scale = 1.0;
                        for (Eigen::Index i = 0; i < count; ++i) {
                            const auto site = static_cast<std::size_t>(i);
                            norm = std::max(norm, wide_math::abs(residuals[site].c[k] * normalization));
                            scale = std::max(scale, wide_math::abs(site_series[site].c[k] * normalization / site_series[site].c[0]));
                            rows(i) = static_cast<double>(-residuals[site].c[k]);
                        }
                        entry.probe("association_refinement_residual", {static_cast<double>(k), static_cast<double>(norm), static_cast<double>(scale)}, used);
                        if (norm <= std::numeric_limits<Base>::epsilon() * scale) break;
                        if (used++ >= static_cast<std::size_t>(policy::ASSOCIATION_MAX_ITERATIONS.value))
                            throw std::runtime_error("extended association residual/implicit-coefficient floor unresolved");
                        coefficient = factorization.solve(rows);
                        bool changed = false;
                        for (Eigen::Index i = 0; i < count; ++i) {
                            auto& x = site_series[static_cast<std::size_t>(i)].c[k];
                            const Base next = x + Base(coefficient(i));
                            if (!wide_math::isfinite(next) || (k == 0 && !(next > 0.0 && next <= 1.0)))
                                throw std::runtime_error("extended association site domain");
                            changed = changed || next != x; x = next;
                        }
                        if (!changed) throw std::runtime_error("extended association residual cannot be refined");
                    }
                }
                // Coupled achieved site defects J0 e_k = -r_k - sum_{l=1}^k J_l e_{k-l}.
                residuals = association_residuals(rho, composition, parameters_, site_series, delta);
                const auto bonded = bonded_densities(composition, parameters_, site_series, delta);
                std::vector<Series> error(sites.size()), jac(sites.size() * sites.size());
                for (std::size_t a = 0; a < sites.size(); ++a) for (std::size_t b = 0; b < sites.size(); ++b) {
                    jac[a * sites.size() + b] = site_series[a] * rho * Base(NUMBER_DENSITY_PER_ANGSTROM3)
                        * composition[parameters_.site_component[b]] * Base(parameters_.site_multiplicity[b]) * delta[a * sites.size() + b];
                    if (a == b) jac[a * sites.size() + b] += Series(1) + rho * Base(NUMBER_DENSITY_PER_ANGSTROM3) * bonded[a];
                }
                for (std::size_t k = 0; k <= N; ++k) {
                    for (std::size_t a = 0; a < sites.size(); ++a) {
                        Base r = residuals[a].c[k];
                        for (std::size_t l = 1; l <= k; ++l) for (std::size_t b = 0; b < sites.size(); ++b)
                            r += jac[a * sites.size() + b].c[l] * error[b].c[k - l];
                        rows(static_cast<Eigen::Index>(a)) = static_cast<double>(-r);
                    }
                    coefficient = factorization.solve(rows);
                    for (Eigen::Index i = 0; i < count; ++i) error[static_cast<std::size_t>(i)].c[k] = Base(coefficient(i));
                    Base defect = 0;
                    for (std::size_t a = 0; a < sites.size(); ++a) {
                        Base r = residuals[a].c[k];
                        for (std::size_t l = 0; l <= k; ++l) for (std::size_t b = 0; b < sites.size(); ++b)
                            r += jac[a * sites.size() + b].c[l] * error[b].c[k - l];
                        defect = std::max(defect, wide_math::abs(r));
                    }
                    entry.probe("association_linear_defect", {double(k), double(defect)});
                }
                for (std::size_t a = 0; a < sites.size(); ++a) {
                    central_sites_[a] = site_series[a].c[0];
                    const auto b = total * composition[parameters_.site_component[a]] * Base(parameters_.site_multiplicity[a]) * residuals[a] / site_series[a];
                    floor += wide_math::abs(b.c[1] * error[a].c[0]) + wide_math::abs(b.c[0] * error[a].c[1]);
                }
            } else {
            // Order k of F(V(t), X(t)) = 0 is F_X X_k + (the residual with X_k = 0): one solve per order.
            for (std::size_t k = 1; k < N; ++k) {
                residuals = association_residuals(rho, composition, parameters_, site_series, delta);
                for (Eigen::Index s = 0; s < count; ++s) rows(s) = -residuals[static_cast<std::size_t>(s)].c[k];
                coefficient = factorization.solve(rows);
                for (Eigen::Index s = 0; s < count; ++s) site_series[static_cast<std::size_t>(s)].c[k] = coefficient(s);
                if (k > 1) continue;
                for (Eigen::Index s = 0; s < count; ++s) {
                    const auto site = static_cast<std::size_t>(s);
                    cross(s) = -total_ * composition_[static_cast<std::size_t>(parameters_.site_component[site])]
                        * parameters_.site_multiplicity[site] / sites[site] * residuals[site].c[1];
                }
                unit = factorization.solve(Eigen::VectorXd::Ones(count));
                floor = wide_math::abs(cross.dot(unit)) * std::max(record.residual, record.floor.value_or(0.0));
            }
            }
            residuals = association_residuals(rho, composition, parameters_, site_series, delta);
        }
        const auto electrostatics = [&] {
            if constexpr (std::is_same_v<Base, double>) return detail::supported_electrostatics(temperature_, rho, composition_, composition_, parameters_);
            else return detail::supported_electrostatics(Series(temperature_), rho, composition, narrow_values(composition_.begin(), composition_.end()), parameters_);
        }();
        const auto contributions = [&] {
            if constexpr (std::is_same_v<Base, double>) return residual_helmholtz(temperature_, rho, composition_, parameters_, sigma_,
                diameter_, site_series, electrostatics.first.debye_huckel, electrostatics.second);
            else return residual_helmholtz(Series(temperature_), rho, composition, parameters_,
                std::vector<Series>(sigma_.begin(), sigma_.end()), std::vector<Series>(diameter_.begin(), diameter_.end()),
                site_series, electrostatics.first.debye_huckel, electrostatics.second);
        }();
        if (sites_out != nullptr) {
            sites_out->clear();
            for (const auto& x : site_series) sites_out->insert(sites_out->end(), x.c.begin(), x.c.begin() + N);
        }
        const auto potential = detail::stationary_association(contributions.association, composition, parameters_, residuals);
        auto result = total * (contributions.hard_chain + contributions.dispersion + contributions.debye_huckel
                         + contributions.born + potential + contributions.dipolar);
        if (!parameters_.cation_mixing.empty()) result += total * contributions.cation_mixing;
        return result;
    }

private:
    const HelmholtzModel& model_;
    std::conditional_t<std::is_same_v<Base, double>, const Parameters&, ParameterRecord<Base>> parameters_;
    Base temperature_, total_ = 0.0;
    bool pairs_ = true;
    mutable std::vector<Base> central_sites_;
    mutable std::size_t used_ = 0;
    std::vector<Base> composition_, sigma_, diameter_, amounts_;
};

PressurePoint at_density(const PressureSeries<>& series, double temperature, double pressure,
                         const std::vector<double>& composition, double density,
                         std::size_t& evaluations, DerivativeOrder order) {
    auto* partition = series.partition ? &*series.partition : nullptr;
    if (partition) if (const auto at = partition->samples.find(density); at != partition->samples.end()) return at->second;
    if (partition && partition->active) {
        if (partition->evaluations == DensityPartition::cap)
            throw std::runtime_error("enumeration evaluation cap exhausted");
        ++partition->evaluations;
    }
    ++evaluations;
    const double rt = GAS_CONSTANT_J_PER_MOL_K * temperature;
    const double scale = density * GAS_CONSTANT_J_PER_MOL_K * temperature;  // z - P/(rho R T), as recorded
    if (order == DerivativeOrder::PressureSlope) {
        double floor = 0.0;
        const auto energy = series.residual<2>(density, floor);
        if (!std::all_of(energy.c.begin(), energy.c.end(), [](double v) { return std::isfinite(v); }))
            throw std::runtime_error("nonfinite density derivatives");
        const double calculated = rt * (density - energy.c[1]);
        const PressurePoint point{calculated / scale - pressure / scale, calculated,
                rt * (1.0 + 2.0 * energy.c[2] / (density * density)), floor};
        if (partition && partition->active) {
            if (!std::isfinite(point.residual) || !std::isfinite(point.slope)) throw std::runtime_error("nonfinite enumeration sample");
            partition->samples.emplace(density, point);
        }
        return point;
    }
    const auto derivatives = series.model().evaluate(temperature, 1.0 / density, composition, order);
    if (!std::all_of(derivatives.gradient.begin(), derivatives.gradient.end(), [](double v) { return std::isfinite(v); })
        || !std::all_of(derivatives.hessian.begin(), derivatives.hessian.end(), [](double v) { return std::isfinite(v); }))
        throw std::runtime_error("nonfinite density derivatives");
    const double calculated = rt * (density - derivatives.gradient[1]);
    return {calculated / scale - pressure / scale, calculated, 0.0, derivatives.gradient_floor};
}

double residual_derivative(
    const PressurePoint& point, double density, double pressure, double temperature
) {
    return (point.slope * density - (point.pressure - pressure))
        / (density * density * GAS_CONSTANT_J_PER_MOL_K * temperature);
}

bool bracket_closed(const PressurePoint& point, double derivative, double width, double tolerance) {
    return point.residual == 0.0 || (std::abs(point.residual) <= tolerance && std::isfinite(derivative)
        && derivative != 0.0 && width * std::abs(derivative) <= tolerance);
}

bool association_failure(const std::exception& error) {
    const auto* failure = dynamic_cast<const core::ConvergenceFailure*>(&error);
    // Who threw, not how it stopped. Enumerating stops silently narrowed this when
    // `classify` gained `stalled`: the scan propagated instead of skipping the point.
    return failure != nullptr
        && failure->record.criterion && std::string_view(failure->record.criterion->name) == policy::ASSOCIATION_RESIDUAL.name;
}

// Safeguarded Newton from the bracket endpoint with smaller absolute pressure error.
// Both endpoint evaluations enclose one crossing in [lower, upper]. The caller
// isolates a mechanically stable interval, or supplies an exact scan-node root.
DensityRoot refine(
    const PressureSeries<>& model, double temperature, double pressure,
    const std::vector<double>& composition, double lower, PressurePoint low,
    double upper, PressurePoint high, core::Scope& parent_scope, std::size_t& evaluations
) {
    const auto maximum = static_cast<std::size_t>(policy::DENSITY_MAX_ITERATIONS.value);
    const double tolerance = model.solver_options().value(policy::density);
    const bool from_lower = std::abs(low.pressure - pressure) < std::abs(high.pressure - pressure);
    double density = from_lower ? lower : upper;
    PressurePoint point = from_lower ? low : high;
    core::ResidualHistory history;
    history.push(std::abs(point.residual));
    std::size_t iterations = 0;
    for (; iterations < maximum; ++iterations) {
        const double derivative = residual_derivative(point, density, pressure, temperature);
        const bool residual_met = std::abs(point.residual) <= tolerance;
        if (bracket_closed(point, derivative, upper - lower, tolerance)) break;
        // Once the residual passes, the bracket must still witness the tolerance/amplification
        // radius: step a quarter radius past the Newton root toward the far end (bisection fallback).
        double candidate = policy::DENSITY_BISECTION_FRACTION.value * (lower + upper);
        if (std::isfinite(derivative) && derivative != 0.0) {
            const double far = upper - density > density - lower ? upper : lower;
            candidate = density - point.residual / derivative
                + (residual_met ? std::copysign(0.25 * tolerance / std::abs(derivative), far - density) : 0.0);
        }
        if (!(candidate > lower && candidate < upper) || !std::isfinite(candidate)) {
            candidate = policy::DENSITY_BISECTION_FRACTION.value * (lower + upper);
        }
        parent_scope.probe(
            "density_newton_step", {density, point.residual, derivative, candidate},
            static_cast<std::uint32_t>(iterations)
        );
        const PressurePoint next = at_density(model, temperature, pressure, composition, candidate,
                                              evaluations, DerivativeOrder::PressureSlope);
        history.push(std::abs(next.residual));
        if (low.residual * next.residual <= 0.0) {
            upper = candidate;
            high = next;
        } else {
            lower = candidate;
            low = next;
        }
        const bool from_lower = std::abs(low.pressure - pressure) < std::abs(high.pressure - pressure);
        density = from_lower ? lower : upper;
        point = from_lower ? low : high;
    }
    const double derivative = residual_derivative(point, density, pressure, temperature);
    parent_scope.probe("density_root_bracket", {lower, upper, derivative, point.residual});
    core::Convergence record{
        std::abs(point.residual), model.solver_options().criterion(policy::density),
        point.gradient_floor / density,
        std::isfinite(derivative) && derivative != 0.0
            ? std::optional<double>(1.0 / std::abs(derivative)) : std::nullopt,
        history.contraction(), iterations + 1, evaluations, core::Stop::stalled, model.solver_options().owners[policy::density]
    };
    record.stop = core::classify(record.residual, *record.criterion, record.floor,
                                 record.contraction, iterations == maximum,
                                 record.amplification);
    if (iterations == maximum && record.stop == core::Stop::converged)
        throw std::runtime_error("density root bracket unresolved within the search budget");
    if (record.stop == core::Stop::nonfinite
        || (iterations == maximum && record.stop != core::Stop::converged)) {
        throw core::ConvergenceFailure(record);
    }
    const DensityRoot root{density, point.pressure - pressure, point.slope, iterations + 1,
                           std::nullopt, DensityBranch::Unspecified, record};
    parent_scope.probe("density_root", {density, point.residual, point.slope, point.slope > 0.0 ? 1.0 : 0.0},
                       static_cast<std::uint32_t>(iterations));
    return root;
}

// Once the single mechanical loop is exposed, each cell contains at most one fold.
// Isolate its positive-slope interval before refining one strict pressure crossing.
// Exact scan-node roots belong to scan_roots, not to either adjacent cell.
// Pressure rising through the target between nonpositive end slopes hides a positive slope
// (mean value theorem): bisect toward the crossing until a stable sample splits the cell into
// two fold cells, or fail. An exact stable root at that sample is this cell's one root.
void scan_cell(
    const PressureSeries<>& model, double temperature, double pressure,
    const std::vector<double>& composition, double lower, PressurePoint low,
    double upper, PressurePoint high, core::Scope& scope, std::size_t& evaluations,
    std::vector<DensityRoot>& roots, std::vector<std::pair<std::size_t, int>>* pattern = nullptr, std::size_t owner = 0
) {
    std::size_t crossing_position = pattern ? pattern->size() : 0;
    if (low.slope <= 0.0 && high.slope <= 0.0) {
        if (pattern && low.residual < 0.0 && high.residual > 0.0)
            throw std::runtime_error("enumeration negative-slope cell hides an upward crossing");
        if (pattern && low.residual > 0.0 && high.residual < 0.0) pattern->emplace_back(owner, -1);
        if (pattern || !(low.residual < 0.0 && high.residual > 0.0)) return;
        for (std::size_t i = 0;; ++i) {
            const double density = lower + policy::DENSITY_BISECTION_FRACTION.value * (upper - lower);
            if (i == static_cast<std::size_t>(policy::DENSITY_MAX_ITERATIONS.value) || !(density > lower && density < upper))
                throw std::runtime_error("density hidden crossing unresolved within the search budget");
            const auto point = at_density(model, temperature, pressure, composition, density,
                                          evaluations, DerivativeOrder::PressureSlope);
            scope.probe("density_hidden_crossing", {density, point.residual, point.slope});
            if (point.slope > 0.0 && point.residual == 0.0)
                roots.push_back(refine(model, temperature, pressure, composition, density, point,
                                       density, point, scope, evaluations));
            if (point.slope > 0.0 && point.residual != 0.0) {
                scan_cell(model, temperature, pressure, composition, lower, low, density, point, scope, evaluations, roots);
                scan_cell(model, temperature, pressure, composition, density, point, upper, high, scope, evaluations, roots);
            }
            if (point.slope > 0.0) return;
            if (point.residual > 0.0) { upper = density; high = point; }
            else { lower = density; low = point; }
        }
    }
    if ((low.slope > 0.0) != (high.slope > 0.0)) {
        const bool stable_left = low.slope > 0.0;
        double a = lower, b = upper;
        PressurePoint pa = low, pb = high;
        // The stable side's pressure runs monotonically from its cell end toward the fold, so a
        // cell end already past the target leaves that side rootless without locating the fold.
        const double margin = model.solver_options().value(policy::density) * GAS_CONSTANT_J_PER_MOL_K * temperature * upper;
        bool bracketed = !pattern && (stable_left ? low.pressure > pressure + margin : high.pressure < pressure - margin);
        for (std::size_t i = 0; i < static_cast<std::size_t>(policy::DENSITY_MAX_ITERATIONS.value)
                 && !bracketed && pa.slope != 0.0 && pb.slope != 0.0; ++i) {
            const double density = a + policy::DENSITY_BISECTION_FRACTION.value * (b - a);
            if (!(density > a && density < b)) break;
            const auto point = at_density(model, temperature, pressure, composition, density,
                                          evaluations, DerivativeOrder::PressureSlope);
            scope.probe("density_spinodal", {density, point.residual, point.slope});
            if ((pa.slope > 0.0) == (point.slope > 0.0)) { a = density; pa = point; }
            else { b = density; pb = point; }
            // P rises across the stable side: a stable point past the target already
            // brackets the unique root; only a rootless side needs the spinodal itself.
            bracketed = !pattern && point.slope > 0.0 && std::abs(point.residual) > model.solver_options().value(policy::density)
                && (stable_left ? low.residual < 0.0 : high.residual > 0.0) && (point.residual > 0.0) == stable_left;
        }
        if (!bracketed && pa.slope != 0.0 && pb.slope != 0.0 && std::nextafter(a, b) < b)
            throw std::runtime_error("density spinodal unresolved within the search budget");
        const bool use_a = pa.slope == 0.0 || (pb.slope != 0.0 && stable_left);
        const double density = use_a ? a : b;
        const auto point = use_a ? pa : pb;
        if (std::abs(point.residual) <= model.solver_options().value(policy::density))
            throw std::runtime_error("density branch unresolved at spinodal pressure within residual tolerance");
        if (pattern && low.slope != 0.0 && high.slope != 0.0) {
            if (!stable_left && low.residual > 0.0 && point.residual < 0.0) pattern->emplace_back(owner, -1);
            pattern->emplace_back(owner, stable_left ? -2 : 2);
            if (stable_left && point.residual > 0.0 && high.residual < 0.0) pattern->emplace_back(owner, -1);
            if (!stable_left) crossing_position = pattern->size();
        }
        if (stable_left) { upper = density; high = point; }
        else { lower = density; low = point; }
    }
    if (pattern && low.residual > 0.0 && high.residual < 0.0) pattern->emplace_back(owner, -1);
    if (low.residual < 0.0 && high.residual > 0.0) {
        if (pattern) pattern->insert(pattern->begin() + static_cast<std::ptrdiff_t>(crossing_position), {owner, 1});
        scope.probe("density_bracket", {lower, upper, low.residual, high.residual});
        roots.push_back(refine(model, temperature, pressure, composition, lower, low,
            upper, high, scope, evaluations));
    }
}


void search_slope(
    const PressureSeries<>& model, double temperature, double pressure,
    const std::vector<double>& composition, std::vector<ScanPoint>& points,
    core::Scope& scope, std::size_t& evaluations,
    const std::function<PressurePoint(double)>& evaluate = {},
    const std::function<void(double, double, const std::string&)>& unresolved = {}
) {
    std::size_t lowest = 0;
    for (std::size_t i = 1; i + 1 < points.size(); ++i) {
        if (points[i].linked && points[i + 1].linked
            && points[i].point.slope <= std::min(points[i - 1].point.slope, points[i + 1].point.slope)
            && (lowest == 0 || points[i].point.slope < points[lowest].point.slope)) lowest = i;
    }
    if (lowest != 0 && points[lowest].point.slope > 0.0) {
        double a = std::log(points[lowest - 1].density), b = std::log(points[lowest + 1].density);
        double x = std::log(points[lowest].density), slope_x = points[lowest].point.slope;
        double slope_a = points[lowest - 1].point.slope, slope_b = points[lowest + 1].point.slope;
        // A positive lower bound for the interpolating quadratic stops this model-based search.
        // It is not a bound on unsampled EOS curvature, even near criticality in log density.
        try {
            for (std::size_t searched = 0;; ++searched) {
                if (!(b - a > std::numeric_limits<double>::epsilon() * (1.0 + std::abs(x))))
                    throw std::runtime_error("density slope search bracket collapsed before resolution");
                const double curvature = 2.0 * (slope_a / ((a - x) * (a - b)) + slope_x / ((x - a) * (x - b))
                                                + slope_b / ((b - a) * (b - x)));
                if (slope_x - 0.5 * std::max(curvature, 0.0) * (b - a) * (b - a) > 0.0) break;
                if (searched == static_cast<std::size_t>(policy::DENSITY_MAX_ITERATIONS.value))
                    throw std::runtime_error("density slope minimum unresolved within the search budget");
                const double u = x - a > b - x ? x - GOLDEN_SECTION * (x - a) : x + GOLDEN_SECTION * (b - x);
                const ScanPoint trial{std::exp(u), evaluate ? evaluate(std::exp(u)) : at_density(model, temperature, pressure, composition,
                    std::exp(u), evaluations, DerivativeOrder::PressureSlope), true};
                scope.probe("density_slope_search", {trial.density, trial.point.residual, trial.point.slope});
                if (trial.point.slope < 0.0) {
                    points.insert(std::upper_bound(points.begin(), points.end(), trial,
                        [](const ScanPoint& l, const ScanPoint& r) { return l.density < r.density; }), trial);
                    break;
                }
                if (trial.point.slope < slope_x) {
                    (u < x ? slope_b : slope_a) = slope_x;
                    (u < x ? b : a) = x;
                    x = u;
                    slope_x = trial.point.slope;
                } else {
                    (u < x ? slope_a : slope_b) = trial.point.slope;
                    (u < x ? a : b) = u;
                }
            }
        } catch (const std::runtime_error& error) {
            if (unresolved) unresolved(std::exp(a), std::exp(b), error.what());
            throw;
        }
    }
}

std::vector<DensityRoot> partition_roots(
    const PressureSeries<>& model, double temperature, double pressure, const std::vector<double>& composition,
    const std::vector<ScanPoint>& points, core::Scope& scope, std::size_t& evaluations,
    std::vector<std::pair<std::size_t, int>>* events = nullptr
) {
    std::vector<DensityRoot> roots;
    for (std::size_t i = 0; i < points.size(); ++i) {
        const auto& node = points[i];
        if (i != 0 && points[i].linked) {
            scan_cell(model, temperature, pressure, composition, points[i - 1].density,
                      points[i - 1].point, points[i].density, points[i].point, scope, evaluations,
                      roots, events, 2 * i - 1);
        }
        if (events && node.point.residual == 0.0) events->emplace_back(2 * i, node.point.slope > 0.0 ? 1 : -1);
        if (events && node.point.slope == 0.0) {
            if (std::abs(node.point.residual) <= model.solver_options().value(policy::density))
                throw std::runtime_error("enumeration target unresolved at an exact spinodal");
            events->emplace_back(2 * i, 4);
        }
        if (node.point.residual == 0.0 && node.point.slope > 0.0)
            roots.push_back(refine(model, temperature, pressure, composition, node.density,
                node.point, node.density, node.point, scope, evaluations));

    }
    const auto witness = std::find_if(points.begin(), points.end(), [](const auto& node) {
        return std::isfinite(node.point.slope) && node.point.slope < 0.0;
    });
    if (witness != points.end())
        for (auto& root : roots) root.negative_slope_density = witness->density;
    return roots;
}

std::vector<double> log_nodes(double log_lower, double log_upper, std::size_t intervals) {
    std::vector<double> nodes;
    for (std::size_t i = 0; i <= intervals; ++i)
        nodes.push_back(std::exp(log_lower + static_cast<double>(i) / intervals * (log_upper - log_lower)));
    return nodes;
}

// Deterministic log-density grid, then the cells. A loop with positive slopes at both ends
// of its cell is invisible to the cells; with one mechanical loop P_rho has one interior
// minimum, so a positive grid minimum of the slope is bracketed and searched by golden
// section in log density until a negative slope exposes the loop or the quadratic slope model
// estimates its minimum positive; this finite search is not a completeness certificate.
std::vector<DensityRoot> scan_roots(
    const PressureSeries<>& model, double temperature, double pressure,
    const std::vector<double>& composition, double log_lower, double log_upper,
    std::size_t intervals, core::Scope& scope, std::size_t& evaluations
) try {
    std::vector<DensityRoot> local;
    if (model.partition) {
        auto& part = *model.partition;
        part.lower = std::max(std::exp(log_lower), part.lower); part.upper = std::min(std::exp(log_upper), part.upper);
        part.active = true;
        if (part.lower < part.upper) {
            auto nodes = log_nodes(std::log(part.lower), std::log(part.upper), static_cast<std::size_t>(policy::DENSITY_BRACKET_INTERVALS.value));
            nodes.front() = part.lower; nodes.back() = part.upper;
            for (const auto& [rho, point] : part.samples) if (rho > part.lower && rho < part.upper) nodes.push_back(rho);
            std::sort(nodes.begin(), nodes.end()); nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
            std::vector<std::pair<std::size_t, int>> previous;
            for (bool first = true;; first = false) {
                const auto missing = std::count_if(nodes.begin(), nodes.end(), [&](double rho) {
                    return part.samples.find(rho) == part.samples.end(); });
                if (part.evaluations + missing > DensityPartition::cap)
                    throw std::runtime_error("enumeration finer mesh cannot fit within evaluation cap");
                std::vector<ScanPoint> points;
                for (double rho : nodes) points.push_back({rho, at_density(model, temperature, pressure, composition, rho, evaluations, DerivativeOrder::PressureSlope), true});
                std::vector<std::pair<std::size_t, int>> events;
                part.intervals = nodes.size() - 1;
                local = partition_roots(model, temperature, pressure, composition, points, scope, evaluations, &events);
                part.stable = local.size();
                if (std::any_of(local.begin(), local.end(), [](const auto& root) { return !(root.pressure_density_derivative > 0.0); }))
                    throw std::runtime_error("enumeration crossing has nonpositive pressure slope");
                auto parents = events;
                for (auto& event : parents) event.first = event.first % 2 == 0 ? event.first / 2 : 2 * (event.first / 4) + 1;
                scope.probe("density_enumeration", {part.lower, part.upper, static_cast<double>(part.intervals), static_cast<double>(local.size()), static_cast<double>(part.evaluations)});
                if (!first && parents == previous) break;
                previous = std::move(events);
                for (std::size_t i = nodes.size() - 1; i != 0; --i) {
                    const double mid = std::exp(0.5 * (std::log(nodes[i - 1]) + std::log(nodes[i])));
                    if (!(mid > nodes[i - 1] && mid < nodes[i])) throw std::runtime_error("enumeration mesh reached floating-point adjacency before agreement");
                    nodes.insert(nodes.begin() + static_cast<std::ptrdiff_t>(i), mid);
                }
            }
        }
        part.active = false;
    }
    std::vector<ScanPoint> points;
    std::string interruptions;
    bool linked = false;
    if (scope.diagnostics_on()) scope.probe("density_scan_input",
        {temperature, pressure, std::exp(log_lower), std::exp(log_upper), static_cast<double>(intervals)});
    auto nodes = log_nodes(log_lower, log_upper, intervals);
    if (model.partition && model.partition->lower < model.partition->upper) {
        nodes.insert(nodes.end(), {model.partition->lower, model.partition->upper});
        std::sort(nodes.begin(), nodes.end()); nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
    }
    for (std::size_t index = 0; index < nodes.size(); ++index) {
        const double density = nodes[index];
        if (model.partition && density > model.partition->lower && density < model.partition->upper) continue;
        if (model.partition && density == model.partition->upper) linked = false;
        try {
            points.push_back({density, at_density(model, temperature, pressure, composition,
                density, evaluations, DerivativeOrder::PressureSlope), linked});
        } catch (const std::exception& error) {
            const bool association = association_failure(error);
            if (scope.diagnostics_on()) scope.probe(association ? "density_scan_skip" : "density_scan_failure",
                {density, static_cast<double>(index), 0.0}, static_cast<std::uint32_t>(index), error.what());
            if (!association && interruptions.empty()) throw;
            interruptions += "node " + std::to_string(index) + " rho "
                + core::optional_number(density) + ": " + error.what() + "\n";
            if (!association)
                throw DensityStructureError("density scan interrupted by association: " + interruptions, model.first_numerical_failure);
            linked = false;
            continue;
        }
        scope.probe("density_scan_point", {density, points.back().point.residual,
            points.back().point.slope, linked ? 1.0 : 0.0}, static_cast<std::uint32_t>(index));
        linked = true;
    }
    if (!interruptions.empty())
        throw DensityStructureError("density scan interrupted by association: " + interruptions, model.first_numerical_failure);
    if (!model.partition) search_slope(model, temperature, pressure, composition, points, scope, evaluations);
    auto roots = partition_roots(model, temperature, pressure, composition, points, scope, evaluations);
    if (model.partition) {
        // Only exact-zero ordinary nodes already own a local crossing.
        for (const auto& root : local) if (std::none_of(points.begin(), points.end(), [&](const auto& node) {
            return node.density == root.molar_density && node.point.residual == 0.0;
        })) roots.push_back(root);
        const auto witness = std::find_if(model.partition->samples.begin(), model.partition->samples.end(), [](const auto& sample) { return sample.second.slope < 0.0; });
        if (witness != model.partition->samples.end()) for (auto& root : roots) root.negative_slope_density = witness->first;
    }
    return roots;
} catch (const std::exception& error) {
    if (model.partition && dynamic_cast<const DensityStructureError*>(&error) == nullptr) throw model.partition->refused(error.what(), model.first_numerical_failure);
    throw;
}

// ------------------------------------------------------------- single-loop structure

// P, P_rho and P_rhorho at one density from the third-order series.
struct Sample final {
    double density;
    PressurePoint point;
    double curvature;
    std::vector<double> sites;  // site fractions' coefficients in volume (three per site): a neighbour's site-solve start
    bool normal_conversion;
};

// Unresolved structure or an exhausted step budget: the caller scans instead.
struct StructureBroken final {};

// ALG-DENSITY's single mechanical loop, sampled: P_rho has one interior minimum rho_m, so P is concave below rho_m
// (P_rhorho < 0, P_rho falling) and convex above it (P_rhorho > 0, P_rho rising). Each sample is checked against
// that order; a resolved contradiction retains its local partition; both exits throw StructureBroken. Tangents bound P from above on the concave side, from
// below on the convex side.
class Isotherm final {
public:
    Isotherm(const PressureSeries<>& series, double temperature, double pressure, double lower, double upper,
             std::size_t& evaluations, core::Scope& scope)
        : series_(series), temperature_(temperature), pressure_(pressure), lower_(lower), upper_(upper),
          evaluations_(evaluations), scope_(scope) {}
    // A broken structure is traced with its cause before the caller scans.
    [[noreturn]] void broken(const char* cause, double density) const {
        scope_.probe("structured_broken", {density}, 0, cause);
        throw StructureBroken{};
    }

    const policy::SolverOptions& solver_options() const { return series_.solver_options(); }
    void at(double density) {
        if (!(density > lower_ && density < upper_)) broken("outside the domain", density);
        if (samples_.size() == static_cast<std::size_t>(policy::DENSITY_MAX_ITERATIONS.value)) broken("sample budget", density);
        if (std::any_of(samples_.begin(), samples_.end(), [&](const Sample& s) { return s.density == density; })) broken("repeated sample", density);
        ++evaluations_;
        double floor = 0.0;
        Taylor<3> energy;
        std::vector<double> sites, start;
        const Sample* nearest = nullptr;
        for (const auto& s : samples_)
            if (!s.sites.empty() && (!nearest || std::abs(std::log(s.density / density)) < std::abs(std::log(nearest->density / density))))
                nearest = &s;
        if (nearest) {  // the nearest sample's site fractions continued along volume, kept inside (0, 1]
            const double dv = 1.0 / density - 1.0 / nearest->density;
            for (std::size_t i = 0; i + 2 < nearest->sites.size(); i += 3) {
                const double x = nearest->sites[i] + dv * (nearest->sites[i + 1] + dv * nearest->sites[i + 2]);
                start.push_back(x > 0.0 && x <= 1.0 ? x : nearest->sites[i]);
            }
        }
        try { energy = series_.residual<3>(density, floor, start.empty() ? nullptr : &start, &sites); }
        catch (const std::exception& error) { broken(error.what(), density); }
        const double rt = GAS_CONSTANT_J_PER_MOL_K * temperature_, scale = density * rt, squared = density * density;
        const double calculated = rt * (density - energy.c[1]), avv = 2.0 * energy.c[2], avvv = 6.0 * energy.c[3];
        const double fourth = squared * squared, cubed = squared * density, u = avv / squared, one_plus_u = 1.0 + u;
        const double v = avvv / fourth, twice_avv = 2.0 * avv, w = twice_avv / cubed, v_plus_w = v + w;
        const std::array conversion{temperature_, density, rt, energy.c[2], energy.c[3], avv, avvv, squared, fourth, cubed,
                                    u, one_plus_u, v, twice_avv, w, v_plus_w, rt * one_plus_u, -rt * v_plus_w};
        Sample sample{density, {calculated / scale - pressure_ / scale, calculated, rt * one_plus_u, floor},
                      -rt * v_plus_w, std::move(sites), std::all_of(conversion.begin(), conversion.end(), [](double value) { return std::isnormal(value); })};
        if (!std::isfinite(sample.point.residual) || !std::isfinite(sample.point.slope)
            || !std::isfinite(sample.curvature) || sample.curvature == 0.0 || sample.point.slope == 0.0) broken("degenerate sample", density);
        const auto at = std::upper_bound(samples_.begin(), samples_.end(), density,
                                         [](double d, const Sample& s) { return d < s.density; });
        // Concave samples lie below convex ones; P_rho falls across concave samples and rises across convex ones
        // The spacing filter controls slope inspection; order_failure owns rounding resolution.
        const bool concave = sample.curvature < 0.0;
        const auto apart = [&](const Sample& s) { return std::abs(s.density - density) > policy::DENSITY_SLOPE_SEPARATION.value * density; };
        for (auto it = at; it != samples_.end(); ++it)
            if ((concave && it->curvature < 0.0 && apart(*it) && !(it->point.slope < sample.point.slope))
                || (!concave && (it->curvature < 0.0 || (apart(*it) && !(it->point.slope > sample.point.slope))))) order_failure(sample, *it, density, apart(*it));
        for (auto it = samples_.begin(); it != at; ++it)
            if ((!concave && it->curvature > 0.0 && apart(*it) && !(it->point.slope < sample.point.slope))
                || (concave && (it->curvature > 0.0 || (apart(*it) && !(it->point.slope > sample.point.slope))))) order_failure(*it, sample, density, apart(*it));
        samples_.insert(at, std::move(sample));
    }
    const std::vector<Sample>& samples() const { return samples_; }
    Isotherm empty() const { return {series_, temperature_, pressure_, lower_, upper_, evaluations_, scope_}; }
    double pressure() const { return pressure_; }
    double lower() const { return lower_; }
    double upper() const { return upper_; }
    double temperature() const { return temperature_; }
    // The concave side ends below the lowest convex sample and the convex side starts above the highest concave one.
    double concave_end() const {
        for (const auto& s : samples_) if (s.curvature > 0.0) return s.density;
        return upper_;
    }
    double convex_start() const {
        for (auto it = samples_.rbegin(); it != samples_.rend(); ++it) if (it->curvature < 0.0) return it->density;
        return lower_;
    }
    // The sample with the most negative P_rho: inside the loop, nearest its slope minimum.
    std::optional<double> witness() const {
        const Sample* best = nullptr;
        for (const auto& s : samples_) if (s.point.slope < 0.0 && (!best || s.point.slope < best->point.slope)) best = &s;
        return best ? std::optional<double>(best->density) : std::nullopt;
    }
    // Extremum over [from, to] of the lower (concave side) or upper (convex side) envelope of the side's tangents:
    // a bound on P there, at a density where the bound is attained.
    std::pair<double, double> envelope(bool concave_side, double from, double to) const {
        std::vector<const Sample*> lines;
        for (const auto& s : samples_) if ((s.curvature < 0.0) == concave_side) lines.push_back(&s);
        if (lines.empty()) return {concave_side ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity(), 0.5 * (from + to)};
        const auto value = [&](double rho) {
            double bound = concave_side ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity();
            for (const auto* s : lines) {
                const double tangent = s->point.pressure + s->point.slope * (rho - s->density);
                bound = concave_side ? std::min(bound, tangent) : std::max(bound, tangent);
            }
            return bound;
        };
        std::vector<double> candidates{from, to};
        for (std::size_t i = 0; i < lines.size(); ++i)
            for (std::size_t j = i + 1; j < lines.size(); ++j) {
                const double ds = lines[i]->point.slope - lines[j]->point.slope;
                if (ds == 0.0) continue;
                const double rho = (lines[j]->point.pressure - lines[i]->point.pressure
                    + lines[i]->point.slope * lines[i]->density - lines[j]->point.slope * lines[j]->density) / ds;
                if (rho > from && rho < to) candidates.push_back(rho);
            }
        std::pair<double, double> best{value(candidates.front()), candidates.front()};
        for (const double rho : candidates) {
            const double v = value(rho);
            if (concave_side ? v > best.first : v < best.first) best = {v, rho};
        }
        return best;
    }

private:
    [[noreturn]] void order_failure(const Sample& a, const Sample& b, double density, bool apart) const {
        const double rt = GAS_CONSTANT_J_PER_MOL_K * temperature_, g = policy::DENSITY_STRUCTURE_ROUNDOFF.value;
        const auto radii = [&](const Sample& s) {
            const double q = std::abs(s.point.slope), c = std::abs(s.curvature);
            const double twice_rt = 2.0 * rt, slope_norm = twice_rt + q, sum = rt + q, four_sum = 4.0 * sum;
            const double divided = four_sum / s.density, curvature_norm = c + divided;
            const double eq = g * slope_norm, ec = g * curvature_norm;
            const std::array operands{rt, q, c, twice_rt, slope_norm, sum, four_sum, divided, curvature_norm, eq, ec};
            const bool normal = s.normal_conversion && std::all_of(operands.begin(), operands.end(), [](double value) { return std::isnormal(value); });
            return std::tuple{eq, ec, normal};
        };
        const auto [eq_a, ec_a, normal_a] = radii(a);
        const auto [eq_b, ec_b, normal_b] = radii(b);
        const bool convex_a = a.curvature > ec_a, convex_b = b.curvature > ec_b;
        const bool concave_a = a.curvature < -ec_a, concave_b = b.curvature < -ec_b;
        const double reversal = concave_a && concave_b ? b.point.slope - a.point.slope : a.point.slope - b.point.slope;
        const double radius_sum = eq_a + eq_b;
        const bool curvature_order = convex_a && concave_b;
        const bool slope_order = apart && ((concave_a && concave_b) || (convex_a && convex_b))
            && std::isnormal(reversal) && std::isnormal(radius_sum) && reversal > radius_sum;
        if (!(normal_a && normal_b && (curvature_order || slope_order))) broken("curvature or slope order", density);
        scope_.probe("structured_broken", {density}, 0, "curvature or slope order");
        std::ostringstream message;
        message << std::setprecision(17) << "density structure contradiction: " << (curvature_order ? "curvature order" : "slope order")
                << "; T=" << temperature_ << " K; P=" << pressure_ << " Pa; rho_a=" << a.density << "; rho_b=" << b.density
                << " mol/m3; q_a=" << a.point.slope << "; q_b=" << b.point.slope << " Pa m3/mol; c_a=" << a.curvature
                << "; c_b=" << b.curvature << " Pa m6/mol2; E_q_a=" << eq_a << "; E_q_b=" << eq_b
                << " Pa m3/mol; E_c_a=" << ec_a << "; E_c_b=" << ec_b << " Pa m6/mol2";
        series_.partition = DensityPartition{a.density, b.density, message.str(), {}};
        for (const auto& sample : samples_) series_.partition->samples.emplace(sample.density, sample.point);
        for (const auto& sample : {a, b}) series_.partition->samples.emplace(sample.density, sample.point);
        throw StructureBroken{};
    }

    const PressureSeries<>& series_;
    double temperature_, pressure_, lower_, upper_;
    std::size_t& evaluations_;
    core::Scope& scope_;
    std::vector<Sample> samples_;
};

struct SideRoot final {
    bool exists = false;
    std::optional<Sample> lower, upper;  // a bracket of the side's stable root, when found
};

// The next density toward P_t from a sample. On the convex (liquid) side, the root of the local quadratic
// P + P_rho d + P_rhorho d^2 / 2 in the stable form: in Newton's direction and at most twice its length. From
// below, where the stiff liquid isotherm's curvature dominates, Halley's step falls far shorter. Elsewhere
// Halley's step when its curvature term keeps Newton's direction and at most doubles the step, which is when
// (P - P_t) P_rhorho <= P_rho^2, else Newton's: concavity keeps a vapor Newton iterate below its root.
double step_toward(const Isotherm& iso, const Sample& from) {
    const double excess = from.point.pressure - iso.pressure(), slope = from.point.slope;
    const double discriminant = slope * slope - 2.0 * excess * from.curvature;
    if (from.curvature > 0.0 && slope > 0.0 && discriminant >= 0.0)
        return from.density - 2.0 * excess / (slope + std::sqrt(discriminant));
    return from.density + (excess * from.curvature <= slope * slope
        ? -2.0 * excess * slope / (2.0 * slope * slope - excess * from.curvature) : -excess / slope);
}

// From one side of a root: the step toward P_t, or once the residual passes a quarter radius past the root (as
// refine), so the next sample brackets it. `past` points from the sample toward the root.
double approach(const Isotherm& iso, const Sample& from, double past) {
    const double derivative = residual_derivative(from.point, from.density, iso.pressure(), iso.temperature());
    if (std::abs(from.point.residual) <= iso.solver_options().value(policy::density) && std::isfinite(derivative) && derivative != 0.0)
        return from.density - from.point.residual / derivative
            + std::copysign(0.25 * iso.solver_options().value(policy::density) / std::abs(derivative), past);
    return step_toward(iso, from);
}

// A new density in (lo, hi): x itself when it lies inside and off every sample, else the interval's midpoint.
double inside(const Isotherm& iso, double x, double lo, double hi) {
    const auto taken = [&](double rho) {
        return std::any_of(iso.samples().begin(), iso.samples().end(), [&](const Sample& s) { return s.density == rho; });
    };
    return x > lo && x < hi && !taken(x) ? x : 0.5 * (lo + hi);
}

// Refine's stop rule on a bracket [lower, upper] of one up-crossing: the sample nearer P_t meets the residual
// tolerance and the bracket is within the tolerance radius. Otherwise the next density: the step toward P_t from that
// sample, or once its residual passes a quarter radius past the root toward the other end; bisection when the
// step leaves the bracket.
std::optional<double> converge(const Isotherm& iso, const Sample& lower, const Sample& upper) {
    const bool from_lower = std::abs(lower.point.pressure - iso.pressure()) < std::abs(upper.point.pressure - iso.pressure());
    const Sample& near = from_lower ? lower : upper;
    const double derivative = residual_derivative(near.point, near.density, iso.pressure(), iso.temperature());
    const double tolerance = iso.solver_options().value(policy::density);
    const bool met = std::abs(near.point.residual) <= tolerance;
    if (bracket_closed(near.point, derivative, upper.density - lower.density, tolerance)) return std::nullopt;
    double x = !(near.point.slope > 0.0) ? 0.5 * (lower.density + upper.density)
        : met && std::isfinite(derivative) && derivative != 0.0
        ? near.density - near.point.residual / derivative
            + std::copysign(0.25 * tolerance / std::abs(derivative), from_lower ? 1.0 : -1.0)
        : step_toward(iso, near);
    if (!(x > lower.density && x < upper.density)) x = 0.5 * (lower.density + upper.density);
    return x;
}

// Main's fold rule (scan_cell): a side whose pressure extremum lies within the residual tolerance of P_t is
// unresolved, neither rooted nor rootless. A deciding sample is checked by its quadratic estimate of that extremum,
// P - P_t - P_rho^2 / (2 P_rhorho) (a loop sample by P - P_t): within the tolerance, the request goes to the scan,
// which applies the rule at its located spinodal. The tolerance in pressure at a density.
double tolerance_pressure(const Isotherm& iso, double density) {
    return iso.solver_options().value(policy::density) * density * GAS_CONSTANT_J_PER_MOL_K * iso.temperature();
}

const Sample& decided(const Isotherm& iso, const Sample& s) {
    const double reach = s.point.pressure - iso.pressure()
        - (s.point.slope > 0.0 ? s.point.slope * s.point.slope / (2.0 * s.curvature) : 0.0);
    if (std::abs(reach) <= tolerance_pressure(iso, s.density)) iso.broken("pressure extremum within the tolerance", s.density);
    return s;
}

// The concave side holds at most one stable root, the up-crossing below its pressure maximum, which is the
// smallest stable root. A sample at or above P_t on this side (concave, or inside the loop) brackets it with a
// rising concave sample below P_t; the side is rootless when the tangent envelope stays below P_t up to its end.
// To establish only whether the root exists, samples go to the envelope's peak; to locate it, steps toward P_t
// climb from below and converge in the bracket.
SideRoot vapor_side(Isotherm& iso, double start, bool locate) {
    for (;;) {
        const auto& samples = iso.samples();
        const Sample* above = nullptr;
        for (const auto& s : samples)
            if (s.point.residual >= 0.0 && (s.curvature < 0.0 || s.point.slope < 0.0)) { above = &s; break; }
        const Sample* below = nullptr;
        for (const auto& s : samples)
            if (s.point.residual < 0.0 && s.curvature < 0.0 && s.point.slope > 0.0 && (!above || s.density < above->density)) below = &s;
        if (above && above->point.residual == 0.0 && above->point.slope > 0.0) return {true, decided(iso, *above), *above};
        if (above && !locate) return {true, std::nullopt, decided(iso, *above)};
        double x = start;
        if (above && below) {
            const auto next = converge(iso, *below, *above);
            if (!next) return {true, *below, decided(iso, *above)};
            x = *next;
        } else if (above) {
            x = above->point.slope > 0.0 && above->curvature < 0.0 ? approach(iso, *above, -1.0) : 0.5 * above->density;
        } else if (below || std::any_of(samples.begin(), samples.end(), [&](const Sample& s) { return s.density <= start; })) {
            const double end = iso.concave_end();
            const auto [bound, peak] = iso.envelope(true, iso.lower(), end);
            if (bound < iso.pressure() - tolerance_pressure(iso, peak)) return {};
            if (!below) x = 0.25 * samples.front().density;
            else x = inside(iso, locate ? std::min(approach(iso, *below, 1.0), peak)
                                        : below->density - below->point.slope / below->curvature, below->density, end);
        }
        iso.at(x);
    }
}

// The convex side mirrored: its stable root is the up-crossing above its pressure minimum, the largest stable root.
// Without full location, retain a stable convex bracket as existence evidence, never a published root density.
SideRoot liquid_side(Isotherm& iso, double start, bool locate) {
    for (;;) {
        const auto& samples = iso.samples();
        const Sample* below = nullptr;
        for (auto it = samples.rbegin(); it != samples.rend(); ++it)
            if (it->point.residual <= 0.0 && (it->curvature > 0.0 || it->point.slope < 0.0)) { below = &*it; break; }
        const Sample* above = nullptr;
        for (auto it = samples.rbegin(); it != samples.rend(); ++it)
            if (it->point.residual > 0.0 && it->curvature > 0.0 && it->point.slope > 0.0 && (!below || it->density > below->density)) above = &*it;
        if (below && below->point.residual == 0.0 && below->point.slope > 0.0) return {true, decided(iso, *below), *below};
        double x = start;
        if (below && above) {
            if (!locate && below->point.slope > 0.0 && below->curvature > 0.0)
                return {true, decided(iso, *below), *above};
            const auto next = converge(iso, *below, *above);
            if (!next) return {true, decided(iso, *below), *above};
            x = *next;
        } else if (below) {
            x = below->point.slope > 0.0 && below->curvature > 0.0 ? approach(iso, *below, 1.0) : 0.5 * (below->density + iso.upper());
        } else if (above || std::any_of(samples.begin(), samples.end(), [&](const Sample& s) { return s.density >= start; })) {
            const double begin = iso.convex_start();
            const auto [bound, trough] = iso.envelope(false, begin, iso.upper());
            if (bound > iso.pressure() + tolerance_pressure(iso, trough)) return {};
            if (!above) x = 0.5 * (samples.back().density + iso.upper());
            // P_rho grows steeply above the liquid spinodal, so its Newton step creeps: no further than the log midpoint.
            else x = inside(iso, std::max(approach(iso, *above, -1.0), trough), begin, above->density);
        }
        iso.at(x);
    }
}

// A negative slope witnesses the loop. Tangents to P_rho propose a sample, but their positive
// intersection cannot exclude a hidden negative dip under the single-minimum assumption.
void loop_search(Isotherm& iso) {
    while (!iso.witness()) {
        const auto& samples = iso.samples();
        const Sample *a = nullptr, *b = nullptr;
        for (const auto& s : samples) {
            if (s.curvature < 0.0) a = &s;
            else if (!b) b = &s;
        }
        if (!a || !b) iso.broken("loop search without both sides", 0.0);
        const double z = (b->point.slope - a->point.slope + a->curvature * a->density - b->curvature * b->density)
            / (a->curvature - b->curvature);
        const double bound = a->point.slope + a->curvature * (z - a->density);
        if (bound > 0.0) iso.broken("loop absence unresolved", z);
        iso.at(z > a->density && z < b->density ? z : 0.5 * (a->density + b->density));
    }
}

// The scan's witness for a sole root: the first node of scan_roots' log grid with a negative slope. Under the single
// loop those nodes lie between the spinodals, so the first one follows from a binary search over the nodes between
// the highest rising sample below the loop and the lowest loop sample (or the first node above it). Throws
// StructureBroken when no node falls inside the loop: the scan's slope search owns that case.
double grid_witness(Isotherm& iso) {
    const double log_lower = std::log(iso.lower()), log_upper = std::log(iso.upper());
    const auto intervals = static_cast<std::size_t>(policy::DENSITY_BRACKET_INTERVALS.value);
    const auto node = [&](std::size_t index) {
        return std::exp(log_lower + static_cast<double>(index) / intervals * (log_upper - log_lower)); };
    const auto slope_at = [&](double rho) {
        for (const auto& s : iso.samples()) if (s.density == rho) return s.point.slope;
        iso.at(rho);
        for (const auto& s : iso.samples()) if (s.density == rho) return s.point.slope;
        return 0.0;
    };
    const auto loop = std::find_if(iso.samples().begin(), iso.samples().end(), [](const Sample& s) { return s.point.slope < 0.0; });
    if (loop == iso.samples().end()) iso.broken("no sample inside the loop", 0.0);
    const double inside = loop->density;
    double rising = iso.lower();
    for (const auto& s : iso.samples()) if (s.density < inside && s.point.slope > 0.0) rising = s.density;
    std::size_t first = 0, last = 0;  // the nodes in (rising, inside]
    while (first <= intervals && !(node(first) > rising)) ++first;
    last = first;
    while (last <= intervals && node(last) <= inside) ++last;
    // Binary search for the first negative slope among nodes [first, last); `last` (above the lowest loop sample)
    // is the candidate when every one of them is positive.
    std::size_t lo = first, hi = last;
    while (lo < hi) {
        const std::size_t mid = lo + (hi - lo) / 2;
        if (slope_at(node(mid)) < 0.0) hi = mid;
        else lo = mid + 1;
    }
    if (lo > intervals || !(slope_at(node(lo)) < 0.0)) iso.broken("loop between scan nodes", inside);
    return node(lo);
}

// The density at packing fraction eta (the upper bound sits at DENSITY_UPPER_FRACTION of the close-packing density).
double structured_start(const Isotherm& iso, double eta) {
    return eta / policy::CLOSE_PACKING_FRACTION.value * iso.upper() / policy::DENSITY_UPPER_FRACTION.value;
}

// Both sides' stable roots or their absence under the single-loop structure, with the loop witness; nullopt
// when a sample breaks the structure. Brackets are kept for the roots a request needs.
struct Structure final {
    SideRoot vapor, liquid;
    std::optional<double> witness;
};

std::optional<Structure> structured_roots(Isotherm& iso, DensityBranch branch, bool sibling) {
    try {
        const double ideal = iso.pressure() / (GAS_CONSTANT_J_PER_MOL_K * iso.temperature());
        // Starts by packing fraction: dilute below 0.037, the loop's typical slope minimum near 0.133 and a liquid at 0.38.
        const double vapor_start = std::min(ideal, structured_start(iso, 0.037)), liquid_start = structured_start(iso, 0.38);
        const double loop_start = structured_start(iso, 0.133);
        Structure result;
        // Unless a sample already lies inside the loop: a loop sample witnesses the loop and, below P_t, the liquid root.
        const auto loop_sample = [&] {
            if (!iso.witness()) iso.at(loop_start);
        };
        // Every request checks the liquid structure first, including vapor requests.
        const bool locate_liquid = branch != DensityBranch::Vapor || sibling;
        result.liquid = liquid_side(iso, liquid_start, locate_liquid);
        loop_sample();
        const bool locate_vapor = branch != DensityBranch::Liquid || sibling;
        result.vapor = vapor_side(iso, locate_vapor ? vapor_start
            : std::min(2.0 * ideal, structured_start(iso, 0.037)), locate_vapor);
        if (!result.vapor.exists && result.liquid.exists && !locate_liquid)
            result.liquid = liquid_side(iso, liquid_start, true);
        if (result.vapor.exists && !result.liquid.exists) {
            Isotherm fresh = iso.empty();
            result.vapor = vapor_side(fresh, vapor_start, true);
        }
        if (!result.vapor.exists && !result.liquid.exists) iso.broken("neither side has a root", 0.0);
        if (result.vapor.exists && result.liquid.exists && !(result.vapor.upper->density < result.liquid.lower->density))
            iso.broken("side roots out of order", result.vapor.upper->density);
        if (!(result.vapor.exists && result.liquid.exists)) loop_search(iso);
        result.witness = iso.witness();
        if (result.witness && !(result.vapor.exists && result.liquid.exists)) result.witness = grid_witness(iso);
        return result;
    } catch (const StructureBroken&) {
        return std::nullopt;
    }
}

DensityRoot anchor(
    const PressureSeries<>& model, double temperature, double pressure,
    const std::vector<double>& composition, DensityBranch branch, double lower, double upper,
    double density, std::size_t& evaluations
) try {
    core::Scope anchor_scope("density_anchor");
    const double window = policy::DENSITY_ANCHOR_LOG_WINDOW.value, log_anchor = std::log(density);
    const double local_lower = std::exp(std::max(std::log(lower), log_anchor - window));
    const double local_upper = std::exp(std::min(std::log(upper), log_anchor + window));
    const auto maximum = static_cast<std::size_t>(policy::DENSITY_MAX_ITERATIONS.value);
    Isotherm iso(model, temperature, pressure, lower, upper, evaluations, anchor_scope);
    bool certificate_available = true;
    const auto sample = [&](double rho) {
        iso.at(rho);
        return *std::find_if(iso.samples().begin(), iso.samples().end(), [&](const Sample& s) { return s.density == rho; });
    };
    core::ResidualHistory history;
    core::Convergence anchor_record;
    std::size_t iteration = 0;
    for (; iteration < maximum; ++iteration) {
        std::optional<Sample> current;
        if (certificate_available) {
            try { current = sample(density); }
            catch (const StructureBroken&) { certificate_available = false; }
        }
        const auto point = current ? current->point
            : at_density(model, temperature, pressure, composition, density,
                         evaluations, DerivativeOrder::PressureSlope);
        history.push(std::abs(point.residual));
        const double derivative = residual_derivative(point, density, pressure, temperature);
        const double candidate = std::isfinite(derivative) && derivative != 0.0
            ? density - point.residual / derivative
            : density;
        anchor_scope.probe("density_newton_step", {density, point.residual, derivative, candidate},
                           static_cast<std::uint32_t>(iteration));
        anchor_record = {
            std::abs(point.residual), model.solver_options().criterion(policy::density),
            point.gradient_floor / density,
            std::isfinite(derivative) && derivative != 0.0
                ? std::optional<double>(1.0 / std::abs(derivative)) : std::nullopt,
            history.contraction(), iteration + 1, evaluations, core::Stop::stalled, model.solver_options().owners[policy::density]
        };
        if (model.partition || (anchor_record.residual <= anchor_record.criterion->value
            && point.slope > 0.0)) {
            // A positive edge slope certifies the local window only under the single slope minimum.
            // Every unproved case uses the inherited local scan, without restarting Newton.
            const auto certificate = [&]() -> std::optional<DensityRoot> {
                if (model.partition || !current || !std::isfinite(derivative) || derivative == 0.0) return std::nullopt;
                try {
                    const Sample edge = sample(current->curvature > 0.0 ? local_lower : local_upper);
                    if ((edge.curvature > 0.0) != (current->curvature > 0.0) || !(edge.point.slope > 0.0)) return std::nullopt;
                    Sample low = *current, high = *current;
                    if (point.residual != 0.0) {
                        const double past = point.residual < 0.0 ? 1.0 : -1.0;
                        const double endpoint = candidate + past * 0.25 * anchor_record.criterion->value / std::abs(derivative);
                        if (!(endpoint >= local_lower && endpoint <= local_upper)) return std::nullopt;
                        const Sample q = sample(endpoint);
                        if (!(q.point.residual * point.residual < 0.0)) return std::nullopt;
                        (past > 0.0 ? high : low) = q;
                    }
                    return refine(model, temperature, pressure, composition, low.density, low.point, high.density,
                                  high.point, anchor_scope, evaluations);
                } catch (const StructureBroken&) {
                    return std::nullopt;
                }
            };
            std::optional<DensityRoot> selected = certificate();
            if (!selected) {
                auto roots = scan_roots(model, temperature, pressure, composition,
                    std::log(local_lower), std::log(local_upper),
                    static_cast<std::size_t>(policy::DENSITY_ANCHOR_SCAN_INTERVALS.value),
                    anchor_scope, evaluations);
                roots.erase(std::remove_if(roots.begin(), roots.end(),
                    [](const auto& root) { return !(root.pressure_density_derivative > 0.0); }), roots.end());
                anchor_record.evaluations = evaluations;
                if (roots.size() == 1) selected = roots.front();

            }
            if (selected) {
                auto result = *selected;
                result.branch = branch;
                result.convergence.iterations += iteration + 1;
                result.convergence.evaluations = evaluations;
                result.iterations = result.convergence.iterations;
                anchor_scope.finish(result.convergence);
                return result;
            }
            anchor_record.stop = core::Stop::infeasible;
            break;
        }
        if (!std::isfinite(derivative) || derivative == 0.0) {
            // classify reports nonfinite for any nonfinite residual, floor, contraction or
            // amplification; a finite record with a dead derivative is singular.
            anchor_record.stop = core::classify(anchor_record.residual, *anchor_record.criterion,
                anchor_record.floor, anchor_record.contraction, false, anchor_record.amplification)
                == core::Stop::nonfinite ? core::Stop::nonfinite : core::Stop::singular;
            break;
        }
        if (!(candidate > local_lower && candidate < local_upper)
            || !std::isfinite(candidate)) {
            anchor_record.stop = core::Stop::left_domain;
            break;
        }
        density = candidate;
    }
    if (iteration == maximum) {
        anchor_record.stop = core::classify(anchor_record.residual, *anchor_record.criterion,
            anchor_record.floor, anchor_record.contraction, true, anchor_record.amplification);
    }
    anchor_record.cause = "density anchor branch " + std::to_string(static_cast<int>(branch))
        + "; anchor " + core::optional_number(std::exp(log_anchor)) + "; window [" + core::optional_number(local_lower)
        + "," + core::optional_number(local_upper) + "] mol/m3";
    anchor_scope.finish(anchor_record);
    throw core::ConvergenceFailure(anchor_record, model.first_numerical_failure);
} catch (const std::exception& error) {
    if (model.partition && dynamic_cast<const DensityStructureError*>(&error) == nullptr) throw model.partition->refused(error.what(), model.first_numerical_failure);
    throw;
}

}  // namespace

DensityDiagnostics HelmholtzModel::density_diagnostics(
    double temperature, double pressure, const std::vector<double>& composition,
    const std::vector<double>& nodes, double density_tolerance, double log_tolerance,
    std::size_t max_evaluations, double max_seconds, bool adaptive,
    std::optional<double> finite_anchor, double finite_upper, DerivativeOrder evaluation_order
) const {
    core::Scope scope("density_diagnostics");
    validate_composition(composition);
    if (!std::isfinite(temperature) || temperature <= 0.0 || !std::isfinite(pressure) || pressure <= 0.0
        || !std::isfinite(density_tolerance) || density_tolerance <= 0.0
        || !std::isfinite(log_tolerance) || log_tolerance <= 0.0
        || !std::isfinite(max_seconds) || max_seconds <= 0.0 || max_evaluations == 0)
        throw std::invalid_argument("invalid density diagnostic inputs or budget");
    DensityDiagnostics result;
    result.solver_options = solver_options();
    result.evaluation_order = evaluation_order;
    result.geometric_cap = molar_density_upper_bound(temperature, composition);
    if (nodes.size() < 2) throw std::invalid_argument("density diagnostics require at least two declared nodes");
    for (std::size_t i = 0; i < nodes.size(); ++i)
        if (!std::isfinite(nodes[i]) || nodes[i] <= 0.0 || nodes[i] >= result.geometric_cap
            || (i && nodes[i] <= nodes[i - 1])) throw std::invalid_argument("invalid declared density grid");
    if (finite_anchor && (!std::isfinite(*finite_anchor) || *finite_anchor <= 0.0
        || *finite_anchor >= result.geometric_cap || !std::isfinite(finite_upper) || finite_upper <= *finite_anchor))
        throw std::invalid_argument("invalid finite-domain interval");
    if (evaluation_order != DerivativeOrder::Gradient && evaluation_order != DerivativeOrder::PressureSlope
        && evaluation_order != DerivativeOrder::Hessian)
        throw std::invalid_argument("finite-domain evaluation requires Gradient, PressureSlope or Hessian");
    const auto started = std::chrono::steady_clock::now();
    const PressureSeries<> series(*this, temperature, composition);
    auto evaluate = [&](double rho, DerivativeOrder order = DerivativeOrder::PressureSlope) {
        if (result.evaluations >= max_evaluations
            || std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count() >= max_seconds)
            throw DiagnosticBudgetExhausted();
        if (rho >= result.geometric_cap) throw std::runtime_error("geometric domain exit");
        auto point = at_density(series, temperature, pressure, composition, rho, result.evaluations, order);
        if (!std::isfinite(point.pressure) || !std::isfinite(point.slope) || !std::isfinite(point.gradient_floor))
            throw std::runtime_error("nonfinite density evaluation");
        return point;
    };
    auto row = [&](double rho, const PressurePoint& point, double lo, double hi, bool spinodal) {
        DensityDiagnosticPoint value;
        value.molar_density = rho;
        value.packing_fraction = rho * packing_fraction_coefficient(temperature, composition);
        value.pressure = point.pressure;
        value.pressure_density_derivative = point.slope;
        value.mechanical_sign = spinodal ? 0 : (point.slope > 0.0) - (point.slope < 0.0);
        value.lower = lo;
        value.upper = hi;
        value.residual = spinodal ? point.slope / (GAS_CONSTANT_J_PER_MOL_K * temperature) : point.pressure - pressure;
        return value;
    };
    auto refine_crossing = [&](double lo, PressurePoint low, double hi, PressurePoint high, bool spinodal) {
        auto value = [&](const PressurePoint& p) { return spinodal ? p.slope : p.pressure - pressure; };
        std::string stop = "bracket_tolerance";
        try {
            while (hi - lo > density_tolerance || std::log(hi / lo) > log_tolerance) {
                if (value(low) == 0.0 || value(high) == 0.0) break;
                double rho = lo + 0.5 * (hi - lo);
                if (!(rho > lo && rho < hi))
                    throw std::runtime_error("density diagnostic bracket unresolved at adjacent densities");
                const auto point = evaluate(rho);
                if ((value(point) > 0.0) == (value(low) > 0.0)) { lo = rho; low = point; }
                else { hi = rho; high = point; }
            }
        } catch (const std::runtime_error& error) {
            auto record = row(lo, low, lo, hi, spinodal);
            record.stop = "unresolved_bracket";
            record.cause = error.what();
            result.unresolved.push_back(record);
            throw;
        }
        const bool use_low = std::abs(value(low)) <= std::abs(value(high));
        if (value(use_low ? low : high) == 0.0) stop = "exact_numerical_zero";
        auto record = row(use_low ? lo : hi, use_low ? low : high, lo, hi, spinodal);
        record.stop = stop;
        return record;
    };
    std::vector<ScanPoint> points;
    bool linked = false;
    try {
        for (double rho : nodes) {
            try {
                const auto point = evaluate(rho);
                auto sample = row(rho, point, rho, rho, false);
                sample.linked = linked;
                result.samples.push_back(sample);
                if (point.pressure == pressure) {
                    sample.stop = "exact_numerical_zero";
                    result.pressure_roots.push_back(sample);
                }
                if (point.slope == 0.0) {
                    auto fold = row(rho, point, rho, rho, true);
                    fold.stop = "exact_numerical_zero";
                    result.spinodals.push_back(fold);
                }
                if (linked && (points.back().point.pressure - pressure) * (point.pressure - pressure) < 0.0)
                    result.declared_brackets.push_back({points.back().density, rho,
                        points.back().point.pressure - pressure, point.pressure - pressure});
                result.declared_crossings = result.declared_brackets.size();
                points.push_back({rho, point, linked});
                linked = true;
            } catch (const DiagnosticBudgetExhausted&) { throw;
            } catch (const std::runtime_error& error) {
                DensityDiagnosticPoint sample;
                sample.molar_density = sample.lower = sample.upper = rho;
                sample.cause = error.what();
                result.samples.push_back(sample);
                linked = false;
            }
        }
        auto crossings = [&](const std::vector<ScanPoint>& cells) {
            for (std::size_t i = 0; i < cells.size(); ++i) {
                const auto& node = cells[i];
                if (node.point.pressure == pressure && std::none_of(result.pressure_roots.begin(),
                    result.pressure_roots.end(), [&](const auto& root) { return root.molar_density == node.density; })) {
                    auto record = row(node.density, node.point, node.density, node.density, false);
                    record.stop = "exact_numerical_zero";
                    result.pressure_roots.push_back(record);
                }
                if (i && node.linked && (cells[i - 1].point.pressure - pressure) * (node.point.pressure - pressure) < 0.0) {
                    const auto& low = cells[i - 1];
                    // A retained sign bracket or exact numerical zero owns its crossing;
                    // enclosure in this connected partition avoids re-refinement, not distance merging.
                    const bool owned = std::any_of(result.pressure_roots.begin(), result.pressure_roots.end(), [&](const auto& root) {
                        return (root.lower >= low.density && root.upper <= node.density)
                            || (root.stop == "exact_numerical_zero" && root.molar_density > low.density && root.molar_density < node.density);
                    });
                    if (!owned) result.pressure_roots.push_back(refine_crossing(low.density, low.point, node.density, node.point, false));
                }
            }
        };
        crossings(points);
        if (adaptive) {
            search_slope(series, temperature, pressure, composition, points, scope, result.evaluations,
                         [&](double rho) {
                auto point = evaluate(rho);
                auto sample = row(rho, point, rho, rho, false);
                sample.declared = false;
                sample.linked = true;
                result.samples.push_back(sample);
                return point;
            }, [&](double lo, double hi, const std::string& cause) {
                DensityDiagnosticPoint bracket;
                bracket.lower = lo; bracket.upper = hi;
                bracket.stop = "unresolved_slope_search"; bracket.cause = cause;
                result.unresolved.push_back(bracket);
            });
        }
        std::vector<ScanPoint> partitions;
        for (std::size_t i = 0; i < points.size(); ++i) {
            const auto& node = points[i];
            if (i && node.linked && points[i - 1].point.slope * node.point.slope < 0.0) {
                const auto& low = points[i - 1];
                const auto fold = refine_crossing(low.density, low.point, node.density, node.point, true);
                result.spinodals.push_back(fold);
                if (adaptive) partitions.push_back({fold.molar_density, evaluate(fold.molar_density), true});
            }
            partitions.push_back(node);
        }
        if (adaptive) crossings(partitions);
        result.stop = result.pressure_roots.empty() ? "no_sampled_crossing" : "resolved_crossings";
        for (const auto& fold : result.spinodals)
            if (std::abs(fold.pressure - pressure) <= solver_options().value(policy::density)
                * fold.molar_density * GAS_CONSTANT_J_PER_MOL_K * temperature)
                result.stop = "unresolved_tangent_pressure_root";
        if (finite_anchor) {
            double lo = *finite_anchor, hi = lo;
            evaluate(lo, evaluation_order);
            result.finite_domain = DensityDiagnosticPoint{};
            result.finite_domain->lower = result.finite_domain->upper = lo;
            const double cap = std::min(finite_upper, result.geometric_cap);
            bool unavailable = false;
            while (hi < cap) {
                hi = std::min(2.0 * lo, cap);
                result.finite_domain->upper = hi;
                try { evaluate(hi, evaluation_order); lo = hi; }
                catch (const DiagnosticBudgetExhausted&) { throw; }
                catch (const std::runtime_error& error) {
                    result.finite_domain->cause = error.what(); unavailable = true; break;
                }
                result.finite_domain->lower = lo;
            }
            while (unavailable && std::nextafter(lo, hi) < hi) {
                const double rho = lo + 0.5 * (hi - lo);
                try { evaluate(rho, evaluation_order); lo = rho; }
                catch (const DiagnosticBudgetExhausted&) { throw; }
                catch (const std::runtime_error& error) {
                    hi = rho; result.finite_domain->cause = error.what();
                }
                result.finite_domain->lower = lo; result.finite_domain->upper = hi;
            }
            result.finite_domain->molar_density = lo;
            result.finite_domain->stop = unavailable ? "finite_to_unavailable" : "declared_interval_finite";
        }
    } catch (const std::runtime_error& error) {
        result.stop = error.what();
        if (result.finite_domain) result.finite_domain->stop = "unresolved_bracket";
    }
    for (auto* roots : {&result.pressure_roots, &result.spinodals})
        std::sort(roots->begin(), roots->end(), [](const auto& a, const auto& b) { return a.molar_density < b.molar_density; });
    core::retain_first(result.first_numerical_failure, series.first_numerical_failure);
    return result;
}

template <typename Base>
PureIsothermalValues<Base> HelmholtzModel::pure_isothermal(Base temperature, Base volume, Base amount) const {
    validate_point(static_cast<double>(temperature), static_cast<double>(volume), {static_cast<double>(amount)});
    if (parameters_.component_count() != 1 || parameters_.ideal_gas || parameters_.charge_number[0] != 0.0)
        throw std::invalid_argument("pure isothermal projection requires one neutral nonideal component");
    PureIsothermalValues<Base> result;
    auto& d = result.derivative;
    std::vector<Base> sites;
    const PressureSeries<Base> series(*this, temperature, {amount});
    const auto energy = series.template residual<2>(1.0 / volume, d.gradient_floor, nullptr, &sites);
    core::retain_first(d.first_numerical_failure, series.first_numerical_failure);
    if (!std::all_of(energy.c.begin(), energy.c.end(), [](Base v) { return wide_math::isfinite(v); }))
        throw std::runtime_error("nonfinite pure isothermal derivatives");
    const std::size_t width = 3 + active_.size();
    d.value = energy.c[0];
    d.gradient.assign(width, 0.0); d.hessian.assign(width * width, 0.0);
    d.gradient[1] = energy.c[1];
    d.gradient[2] = (d.value - volume * d.gradient[1]) / amount;
    d.volume_second_derivative = 2.0 * energy.c[2];
    d.hessian[width + 1] = d.volume_second_derivative;
    d.hessian[width + 2] = d.hessian[2 * width + 1] = -volume * d.volume_second_derivative / amount;
    d.hessian[2 * width + 2] = volume * volume * d.volume_second_derivative / (amount * amount);
    for (std::size_t i = 0; i < sites.size(); i += 2) d.site_fractions.push_back(static_cast<double>(sites[i]));
    const Base rt = GAS_CONSTANT_J_PER_MOL_K * temperature;
    result.pressure = rt * (amount / volume - d.gradient[1]);
    result.chemical_potential_over_rt = wide_math::log(amount / volume) + d.gradient[2];
    result.pressure_density_derivative = rt * (1.0 + volume * volume * d.volume_second_derivative / amount);
    return result;
}

HelmholtzDerivativeValues<WiderScalar> HelmholtzModel::isothermal(WiderScalar temperature, WiderScalar volume,
    const std::vector<WiderScalar>& amounts) const {
    using Base = WiderScalar;
    std::vector<double> declared;
    for (const auto& amount : amounts) declared.push_back(static_cast<double>(amount));
    validate_point(double(temperature), double(volume), declared);
    HelmholtzDerivativeValues<Base> out;
    const PressureSeries<Base> series(*this, temperature, amounts);
    const auto count = amounts.size() + 1, width = count + 1 + active_.size();
    out.gradient.assign(width, 0); out.hessian.assign(width * width, 0); out.gradient_error.assign(width, 0);
    std::vector<Base> diagonal(count), scales(count, 1), direction(count, 0); scales[0] = volume;
    const auto along = [&](std::size_t i, std::size_t j, Base& floor) {
        std::fill(direction.begin(), direction.end(), 0); direction[i] = scales[i];
        if (j < count) direction[j] = scales[j];
        return series.residual<2>(1 / volume, floor, nullptr, nullptr, direction);
    };
    for (std::size_t i = 0; i < count; ++i) {
        Base floor = 0; const auto energy = along(i, count, floor);
        out.value = energy.c[0]; out.gradient[i + 1] = energy.c[1] / scales[i];
        out.gradient_error[i + 1] = floor / scales[i];
        out.gradient_floor = std::max(out.gradient_floor, out.gradient_error[i + 1]);
        diagonal[i] = energy.c[2]; out.hessian[(i + 1) * width + i + 1] = 2 * diagonal[i] / (scales[i] * scales[i]);
    }
    for (std::size_t i = 0; i < count; ++i) for (std::size_t j = 0; j < i; ++j) {
        Base floor = 0; const auto energy = along(i, j, floor);
        out.hessian[(i + 1) * width + j + 1] = out.hessian[(j + 1) * width + i + 1]
            = (energy.c[2] - diagonal[i] - diagonal[j]) / (scales[i] * scales[j]);
    }
    out.volume_second_derivative = out.hessian[width + 1];
    core::retain_first(out.first_numerical_failure, series.first_numerical_failure);
    return out;
}

template EPC_API PureIsothermalValues<double> HelmholtzModel::pure_isothermal(double, double, double) const;
template PureIsothermalValues<WiderScalar> HelmholtzModel::pure_isothermal(WiderScalar, WiderScalar, WiderScalar) const;

Taylor<3> properties::density_volume_taylor(const HelmholtzModel& model, double temperature, double density,
    const std::vector<double>& composition, double& gradient_floor) {
    return PressureSeries<>(model, temperature, composition).residual<3>(density, gradient_floor);
}

bool properties::regular_density(double density, double slope, const std::optional<core::Convergence>& record) {
    if (!(density > 0.0) || !std::isfinite(density) || !(slope > 0.0) || !std::isfinite(slope)
        || !record || !record->criterion || !record->floor || !record->amplification) return false;
    const auto& c = *record;
    return std::isfinite(c.residual) && c.residual <= c.criterion->value
        && *c.floor >= 0.0 && std::isfinite(*c.floor)
        && *c.amplification > 0.0 && std::isfinite(*c.amplification)
        && (c.criterion->value + *c.floor) * *c.amplification / density
            <= policy::DENSITY_ANCHOR_REGULARITY.value;
}

DensityBranch properties::density_side(DensityBranch branch, double density, std::optional<double> sample) {
    if (branch != DensityBranch::Unspecified || !sample) return branch;
    return density < *sample ? DensityBranch::Vapor : DensityBranch::Liquid;
}

std::optional<std::string> properties::vapor_side(const State& state) {
    if (density_side(state.density_branch, state.molar_density, state.negative_slope_density) != DensityBranch::Vapor) return std::nullopt;
    return "liquid branch not admitted by the finite search: the sole stable root " + core::optional_number(state.molar_density)
        + " mol/m3 lies below the negative-slope sample " + core::optional_number(state.negative_slope_density) + " mol/m3";
}

DensityRoot solve_density(
    const HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& mole_fractions, DensityBranch branch,
    std::optional<double> anchor_density, DensityRoot* sibling_root
) {
    core::Scope scope("solve_density");
    if (anchor_density.has_value() && sibling_root != nullptr)
        throw std::invalid_argument("density sibling requires an unanchored root search");
    if (!std::isfinite(temperature) || temperature <= 0.0) throw std::invalid_argument("invalid temperature " + core::optional_number(temperature) + " K; need T > 0 K");
    if (!std::isfinite(pressure) || pressure <= 0.0) throw std::invalid_argument("invalid pressure");
    model.validate_composition(mole_fractions);
    // EqID: ideal_gas_density
    const double ideal_density = pressure
        / (GAS_CONSTANT_J_PER_MOL_K * temperature);
    if (model.parameters().ideal_gas) {
        if (!(ideal_density > 0.0) || !std::isfinite(ideal_density) || !std::isfinite(GAS_CONSTANT_J_PER_MOL_K * temperature))
            throw std::domain_error("nonfinite ideal-gas density closure");
        if (branch == DensityBranch::Liquid) throw std::invalid_argument("ideal-gas has only a vapor density root");
        if (anchor_density && (!std::isfinite(*anchor_density) || *anchor_density <= 0.0))
            throw std::invalid_argument("invalid density anchor");
        DensityRoot root;
        root.molar_density = ideal_density;
        root.pressure_density_derivative = GAS_CONSTANT_J_PER_MOL_K * temperature;
        root.stable_root_count = 1;
        root.branch = DensityBranch::Vapor;
        core::stopping_criterion(root.convergence, model.solver_options(), policy::density);
        root.convergence.residual = 0.0;
        root.convergence.stop = core::classify(root.convergence.residual, *root.convergence.criterion, std::nullopt, std::nullopt, false);
        scope.finish(root.convergence);
        if (sibling_root) *sibling_root = root;
        return root;
    }
    const double lower = ideal_density * policy::DENSITY_LOWER_FRACTION.value;
    const double upper = model.molar_density_upper_bound(temperature, mole_fractions)
        * policy::DENSITY_UPPER_FRACTION.value;
    if (!(lower > 0.0) || !(lower < upper)) throw std::runtime_error("density domain is empty");
    std::size_t evaluations = 0;
    DensityRoot result;
    const PressureSeries<> series(model, temperature, mole_fractions);
    if (anchor_density.has_value()) {
        try {
            if (branch == DensityBranch::Unspecified) throw std::invalid_argument("density anchor requires explicit branch");
            if (!std::isfinite(*anchor_density) || *anchor_density <= 0.0) {
                throw std::invalid_argument("invalid density anchor");
            }
            if (!(*anchor_density > lower && *anchor_density < upper)) {
                core::Convergence record{0.0, model.solver_options().criterion(policy::density), std::nullopt,
                    std::nullopt, std::nullopt, 0, evaluations, core::Stop::left_domain, model.solver_options().owners[policy::density]};
                record.cause = "density anchor outside " + std::to_string(static_cast<int>(branch))
                    + " domain; anchor=" + core::optional_number(*anchor_density) + "; lower=" + core::optional_number(lower)
                    + "; upper=" + core::optional_number(upper);
                throw core::ConvergenceFailure(record);
            }
            result = anchor(series, temperature, pressure, mole_fractions, branch, lower, upper,
                            *anchor_density, evaluations);
        } catch (const std::exception& error) {
            if (scope.diagnostics_on()) scope.probe("anchor_failure", {*anchor_density}, 0, error.what());
            throw;
        }
    } else {
        Isotherm iso(series, temperature, pressure, lower, upper, evaluations, scope);
        const auto structure = structured_roots(iso, branch, sibling_root != nullptr);
        std::vector<DensityRoot> stable;
        std::size_t count = 0;
        if (structure) {
            const auto refined = [&](const SideRoot& side) {
                return refine(series, temperature, pressure, mole_fractions, side.lower->density, side.lower->point,
                              side.upper->density, side.upper->point, scope, evaluations);
            };
            const bool both = structure->vapor.exists && structure->liquid.exists;
            count = static_cast<std::size_t>(structure->vapor.exists) + static_cast<std::size_t>(structure->liquid.exists);
            if (structure->vapor.exists && (branch == DensityBranch::Vapor || sibling_root != nullptr || !both))
                stable.push_back(refined(structure->vapor));
            if (structure->liquid.exists && (branch == DensityBranch::Liquid || sibling_root != nullptr || !both))
                stable.push_back(refined(structure->liquid));
            for (auto& root : stable) root.negative_slope_density = structure->witness;
        } else {
            const auto roots = scan_roots(series, temperature, pressure, mole_fractions,
                std::log(lower), std::log(upper),
                static_cast<std::size_t>(policy::DENSITY_BRACKET_INTERVALS.value), scope, evaluations);
            std::copy_if(roots.begin(), roots.end(), std::back_inserter(stable), [](const auto& root) {
                return std::isfinite(root.pressure_density_derivative) && root.pressure_density_derivative > 0.0;
            });
            count = stable.size();
            // A sole root takes its own side's value from no samples, as the structured search returns it.
            if (stable.size() == 1 && !series.partition) {
                for (const bool vapor : {true, false}) {
                    try {
                        Isotherm fresh = iso.empty();
                        const auto side = vapor ? vapor_side(fresh, std::min(pressure / (GAS_CONSTANT_J_PER_MOL_K * temperature),
                                                                             structured_start(fresh, 0.037)), true)
                                                : liquid_side(fresh, structured_start(fresh, 0.38), true);
                        if (!side.exists) continue;
                        auto root = refine(series, temperature, pressure, mole_fractions, side.lower->density, side.lower->point,
                                           side.upper->density, side.upper->point, scope, evaluations);
                        // The same root when the two lie within their tolerance radii, (tolerance + floor) / |dr/drho|.
                        const auto radius = [](const DensityRoot& r) {
                            return (r.convergence.criterion->value + r.convergence.floor.value_or(0.0)) * r.convergence.amplification.value_or(0.0); };
                        if (std::abs(root.molar_density - stable.front().molar_density) <= radius(root) + radius(stable.front())) {
                            root.negative_slope_density = stable.front().negative_slope_density;
                            stable.front() = root;
                        }
                        break;
                    } catch (const StructureBroken&) {
                        if (series.partition) {
                            stable = scan_roots(series, temperature, pressure, mole_fractions, std::log(lower), std::log(upper),
                                static_cast<std::size_t>(policy::DENSITY_BRACKET_INTERVALS.value), scope, evaluations);
                            count = stable.size(); break;
                        }
                    }
                }
            }
        }
        scope.probe("stable_root_count", static_cast<double>(count));
        if (series.partition && (stable.empty() || (count > 1 && branch == DensityBranch::Unspecified)))
            throw series.partition->refused("branch ambiguous or not established by finite search; stable=" + std::to_string(count), series.first_numerical_failure);
        if (count > 1 && branch == DensityBranch::Unspecified)
            throw std::runtime_error("density root is ambiguous (stable count "
                                     + std::to_string(count) + ")");
        if (stable.empty()) throw std::runtime_error("no stable density root (count 0) at T = " + core::optional_number(temperature)
            + " K, P = " + core::optional_number(pressure) + " Pa; check the units (K, Pa)");
        std::sort(stable.begin(), stable.end(), [](const auto& left, const auto& right) {
            return left.molar_density < right.molar_density;
        });
        result = branch == DensityBranch::Vapor ? stable.front() : stable.back();
        result.stable_root_count = count;
        result.branch = count == 1 ? DensityBranch::Unspecified : branch;
        if (sibling_root != nullptr) {
            *sibling_root = branch == DensityBranch::Vapor ? stable.back() : stable.front();
            sibling_root->stable_root_count = count;
            sibling_root->branch = count == 1 ? DensityBranch::Unspecified
                : branch == DensityBranch::Vapor ? DensityBranch::Liquid : DensityBranch::Vapor;
        }
    }
    if (series.partition) {
        result.convergence.evaluations = evaluations;
        if (sibling_root) sibling_root->convergence.evaluations = evaluations;
    }
    scope.probe("density_root", result.molar_density);
    core::retain_first(result.first_numerical_failure, series.first_numerical_failure);
    if (sibling_root) core::retain_first(sibling_root->first_numerical_failure, series.first_numerical_failure);
    scope.finish(result.convergence);
    return result;

}

namespace {
// Bell and Deiters (2023), doi:10.1021/acs.iecr.2c02916; start domain verified through this reduced temperature.
constexpr double SUPERANCILLARY_MAX_THETA = 0.999;
double superancillary_series(std::size_t collection, double x) {
    auto first = superancillary::offsets[collection], last = superancillary::offsets[collection + 1];
    while (first + 1 < last) {
        const auto middle = (first + last) / 2;
        if (x < superancillary::segments[middle][0]) last = middle;
        else first = middle;
    }
    const auto& s = superancillary::segments[first];
    const double u = (2.0 * x - s[0] - s[1]) / (s[1] - s[0]);
    double a = 0.0, b = 0.0;
    for (int j = 16; j > 0; --j) { const double previous = a; a = 2.0 * u * a - b + s[j + 2]; b = previous; }
    return u * a - b + s[2];
}
}
std::optional<std::array<double, 2>> superancillary_densities(const HelmholtzModel& model, double temperature) {
    const auto& p = model.parameters();
    const auto zero = [](const auto& v) { return std::all_of(v.begin(), v.end(), [](auto x) { return x == 0; }); };
    if (p.component_count() != 1 || !zero(p.charge_number) || !zero(p.molecular_ion) || !p.site_component.empty()
        || !zero(p.dipole_moment_debye) || !zero(p.quadrupole_moment_debye_angstrom) || !p.sigma_term_component.empty()) return std::nullopt;
    const double m = p.segment_count[0], reduced = temperature / p.dispersion_energy_over_k[0];
    if (!(m >= 1 && m <= 64)) return std::nullopt;
    const double w = 1.0 / m, critical = superancillary_series(0, w);
    const auto& lower = superancillary::lower_temperature[0];
    const double minimum = std::exp(lower[0]) * std::pow(m, lower[1]) * critical;
    const double theta = (reduced - minimum) / (critical - minimum);
    if (!(theta >= 0.0 && theta <= SUPERANCILLARY_MAX_THETA)) return std::nullopt;
    std::size_t interval = 0;
    while (w > superancillary::nodes[interval][0]) ++interval;
    std::array<double, 2> densities{};
    for (std::size_t phase = 0; phase < densities.size(); ++phase) {
        double numerator = 0.0, denominator = 0.0;
        for (std::size_t j = 0; j < 17; ++j) {
            const double value = superancillary_series(1 + interval * 34 + phase * 17 + j, theta);
            const double delta = w - superancillary::nodes[interval][j];
            if (delta == 0.0) { numerator = value; denominator = 1.0; break; }
            const double weight = (j % 2 == 0 ? 1.0 : -1.0) * (j == 0 || j == 16 ? 0.5 : 1.0) / delta;
            numerator += weight * value; denominator += weight;
        }
        densities[phase] = numerator / denominator / (NUMBER_DENSITY_PER_ANGSTROM3 * std::pow(p.sigma_constant[0], 3));
    }
    return densities;
}

}  // namespace epcsaft
