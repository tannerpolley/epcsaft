#pragma once

#include "core/convergence.hpp"
#include "core/numerical_policy.hpp"
#include "core/parameters.hpp"
#include "core/physical_constants.hpp"
#include "core/scalar.hpp"
#include "core/trace.hpp"
#include "core/linear_error.hpp"
#include "contributions/helmholtz.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace epcsaft {

template <typename Scalar, typename ParameterSet, typename Fixed>
std::vector<Scalar> association_delta(
    const Fixed& temperature,
    const Scalar& molar_density,
    const std::vector<Fixed>& mole_fractions,
    const ParameterSet& parameters,
    const std::vector<Fixed>& sigma,
    const std::vector<Fixed>& diameter,
    const std::array<Fixed, 4>* prepared = nullptr
) {
    const std::size_t sites = parameters.site_component.size();
    const auto zeta = packing_moments(molar_density, mole_fractions, parameters, diameter, prepared);
    // Declared laws, then thermal combining rules over their explicit self edges, give energy and
    // volume at this temperature; without them the stored values are read unchanged.
    const auto laws = temperature_law_values(temperature, parameters);
    const bool thermal = !parameters.temperature_law_slot.empty() || !parameters.association_thermal_rule.empty();
    std::vector<Fixed> law_energy, law_volume;
    if (thermal) {
        law_energy.assign(parameters.association_energy_over_k.begin(), parameters.association_energy_over_k.end());
        law_volume.assign(parameters.association_volume.begin(), parameters.association_volume.end());
        for (std::size_t a = 0; a < sites; ++a)
            for (std::size_t b = a; b < sites; ++b) {
                const std::size_t slot = a * sites + b, mirror = b * sites + a;
                if (auto value = temperature_law_value(temperature, parameters, laws, TemperatureLaw::AssociationEnergyOverK, slot, sites))
                    law_energy[slot] = law_energy[mirror] = *value;
                if (auto value = temperature_law_value(temperature, parameters, laws, TemperatureLaw::AssociationVolume, slot, sites))
                    law_volume[slot] = law_volume[mirror] = *value;
            }
        for (const auto& [edge, reverse, self_1, self_2, correction] : parameters.association_thermal_rule) {
            const auto slot = static_cast<std::size_t>(edge), mirror = static_cast<std::size_t>(reverse);
            const auto first = static_cast<std::size_t>(self_1), second = static_cast<std::size_t>(self_2);
            law_energy[slot] = law_energy[mirror] = Fixed(0.5) * (law_energy[first] + law_energy[second]);
            // EqID: association_volume_correction
            Fixed volume = scalar_sqrt(law_volume[first] * law_volume[second]);
            if (correction >= 0) volume *= Fixed(1.0) - laws[static_cast<std::size_t>(correction)];
            if constexpr (std::is_same_v<Fixed, double>)
                if (!(volume >= 0.0)) throw std::invalid_argument(law_domain_failure("association_volume", slot, sites, temperature));
            law_volume[slot] = law_volume[mirror] = volume;
        }
    }
    std::vector<Scalar> delta(sites * sites, Scalar(0.0));
    for (std::size_t a = 0; a < sites; ++a) {
        const auto i = static_cast<std::size_t>(parameters.site_component[a]);
        for (std::size_t b = 0; b < sites; ++b) {
            const std::size_t pair = a * sites + b;
            if (parameters.association_active[pair] == 0) {
                continue;
            }
            const auto j = static_cast<std::size_t>(parameters.site_component[b]);
            // EqID: association_contact
            const Scalar contact = Scalar(1.0) + hard_sphere_contact_departure(
                diameter[i] * diameter[j] / (diameter[i] + diameter[j]), zeta[2], zeta[3]);
            // Wolbach--Sandler pairs (active == 2): the Gross2002 Eq. 3 size
            // factor times sigma_ij^3 is the geometric-mean diameter cubed.
            const Fixed size_ij = (parameters.association_active[pair] == 2
                    ? scalar_sqrt(sigma[i] * sigma[j]) : Fixed(0.5) * (sigma[i] + sigma[j]))
                * (Fixed(1.0) - pair_correction(temperature, parameters, laws, TemperatureLaw::Lij, i, j));
            // EqID: association_strength
            const Fixed energy = thermal ? law_energy[pair] : Fixed(parameters.association_energy_over_k[pair]);
            const Fixed volume = thermal ? law_volume[pair] : Fixed(parameters.association_volume[pair]);
            delta[pair] = contact
                * (scalar_exp(
                       energy / temperature
                   ) - Fixed(1.0))
                * size_ij * size_ij * size_ij
                * volume;
        }
    }
    return delta;
}

template <typename Scalar, typename ParameterSet, typename Fixed>
std::vector<Scalar> association_residuals(
    const Scalar& molar_density,
    const std::vector<Fixed>& mole_fractions,
    const ParameterSet& parameters,
    const std::vector<Scalar>& site_fractions,
    const std::vector<Scalar>& delta
) {
    return association_products(molar_density, mole_fractions, parameters, site_fractions, delta).second;
}

// EqExclude: numerical iteration mechanics
// At vacuum X=1 and the mass-action Jacobian is exactly the identity.
template <typename Scalar, typename Fixed>
std::vector<Scalar> association_density_coefficients(const Scalar& density,
    const std::vector<Fixed>& composition, const Parameters& parameters,
    const std::vector<Scalar>& delta) {
    std::vector<Scalar> sites(parameters.site_count(), Scalar(1.0));
    for (std::size_t k = 1; k <= 2; ++k) {
        const auto rows = association_residuals(density, composition, parameters, sites, delta);
        for (std::size_t a = 0; a < sites.size(); ++a) sites[a].c[k] = -rows[a].c[k];
    }
    return sites;
}

// EqExclude: numerical iteration mechanics
// dF_a/dX_b = X_a rho_N x_c(b) n_b Delta_ab + kron_ab (1 + rho_N B_a) for
// F_a = X_a (1 + rho_N B_a) - 1.
inline Eigen::MatrixXd mass_action_jacobian(
    const std::vector<double>& mole_fractions, const Parameters& parameters,
    const std::vector<double>& site_fractions, double number_density,
    const std::vector<double>& delta, const std::vector<double>& bonded
) {
    const auto sites = static_cast<Eigen::Index>(site_fractions.size());
    Eigen::MatrixXd jacobian(sites, sites);
    for (Eigen::Index a = 0; a < sites; ++a) {
        for (Eigen::Index b = 0; b < sites; ++b)
            jacobian(a, b) = site_fractions[a] * number_density
                * mole_fractions[static_cast<std::size_t>(parameters.site_component[b])]
                * parameters.site_multiplicity[b] * delta[a * sites + b];
        jacobian(a, a) += 1.0 + number_density * bonded[a];
    }
    return jacobian;
}

// Full pivoting: a smallest pivot at or below ASSOCIATION_JACOBIAN_CONDITION times the
// largest marks the Jacobian ill-conditioned. Returns largest/smallest when admitted.
inline std::optional<double> admitted_pivot_ratio(const Eigen::FullPivLU<Eigen::MatrixXd>& factorization) {
    const auto pivots = factorization.matrixLU().diagonal().cwiseAbs();
    const double largest = pivots.maxCoeff(), smallest = pivots.minCoeff();
    if (!std::isfinite(largest) || !(largest > 0.0)
        || smallest <= policy::ASSOCIATION_JACOBIAN_CONDITION.value * largest) return std::nullopt;
    return largest / smallest;
}

// EqExclude: numerical iteration mechanics
// Dense and small: the site count is single digits, so the factorisation never
// shows up next to the tape sweeps. A nearly singular Jacobian returns nothing,
// the signal to fall back to substitution.
inline std::optional<std::vector<double>> solve_linear(
    const Eigen::MatrixXd& matrix,
    const std::vector<double>& right_hand_side,
    double* pivot_ratio = nullptr, std::optional<core::Convergence>* first = nullptr,
    std::optional<core::NumericalRank>* record = nullptr, const std::vector<double>& point = {}, std::size_t iteration = 0
) {
    const Eigen::FullPivLU<Eigen::MatrixXd> factorization(matrix);
    const auto ratio = admitted_pivot_ratio(factorization);
    if (!ratio || factorization.rank() != matrix.cols()) {
        if (record && !*record) {
            auto rank = core::numerical_rank(factorization, "mass_action_site_fractions");
            core::ConditioningObservation warning{"association", "", "pivot", "rank/pivot", "association_pivot_refusal"};
            warning.iteration = iteration; warning.point_kind = "refused_linearization";
            core::observe_conditioning(rank, std::move(warning), point); *record = std::move(rank);
        }
        if (first) core::retain_first(*first, [&] {
            core::Convergence r;
            const auto pivots = factorization.matrixLU().diagonal().cwiseAbs();
            r.residual = Eigen::Map<const Eigen::VectorXd>(right_hand_side.data(), right_hand_side.size()).cwiseAbs().maxCoeff();
            r.amplification = pivots.maxCoeff() / pivots.minCoeff();
            r.rank = std::make_shared<core::NumericalRank>(core::numerical_rank(factorization, "mass_action_site_fractions"));
            r.cause = "association Jacobian is ill-conditioned; largest " + core::optional_number(pivots.maxCoeff())
                + "; smallest " + core::optional_number(pivots.minCoeff()) + "; relative pivot cutoff "
                + core::optional_number(policy::ASSOCIATION_JACOBIAN_CONDITION.value);
            return r;
        });
        return std::nullopt;
    }
    if (pivot_ratio != nullptr) *pivot_ratio = *ratio;
    const auto size = static_cast<Eigen::Index>(right_hand_side.size());
    const Eigen::VectorXd solution = factorization.solve(
        Eigen::Map<const Eigen::VectorXd>(right_hand_side.data(), size)
    );
    return std::vector<double>(solution.data(), solution.data() + size);
}

// EqExclude: numerical iteration mechanics
// `start` replaces the unbonded start X = 1 (a neighbouring solution, inside (0, 1]); `strength` supplies
// association_delta at this state; `pair_start` takes the analytic pair step from X = 1 whatever the association
// strength. None changes the solved system or its tolerance.
inline std::vector<double> solve_site_fractions(
    double temperature,
    double molar_density,
    const std::vector<double>& mole_fractions,
    const Parameters& parameters,
    const std::vector<double>& sigma,
    const std::vector<double>& diameter,
    core::Convergence* record = nullptr,
    const std::vector<double>* start = nullptr,
    const std::vector<double>* strength = nullptr,
    bool pair_start = false,
    const policy::SolverOptions& options = {}, std::optional<core::Convergence>* first_numerical_failure = nullptr,
    const std::array<double, 4>* prepared = nullptr
) {
    std::optional<core::Convergence> local_first;
    auto& first = first_numerical_failure ? *first_numerical_failure : local_first;
    if (prepared) (void)packing_fraction_from_geometry(molar_density, (*prepared)[3]);
    else (void)packing_fraction_from_geometry(molar_density, mole_fractions, parameters, diameter);
    const std::size_t sites = parameters.site_component.size();
    if (sites == 0) {
        if (record != nullptr) {
            *record = {};
            core::stopping_criterion(*record, options, policy::association);
            record->floor = 0.0;
            record->amplification = 1.0;
            record->stop = core::classify(
                record->residual, *record->criterion, record->floor,
                record->contraction, false, record->amplification
            );
        }
        return {};
    }
    core::Scope scope("association_newton");
    const auto computed = strength != nullptr ? std::vector<double>{}
        : association_delta(temperature, molar_density, mole_fractions, parameters, sigma, diameter, prepared);
    const auto& delta = strength != nullptr ? *strength : computed;
    const double number_density =
        molar_density * NUMBER_DENSITY_PER_ANGSTROM3;
    std::vector<double> fractions = start != nullptr ? *start : std::vector<double>(sites, 1.0);
    core::ResidualHistory history;
    std::size_t evaluations = 0;
    std::optional<double> last_pivot;
    std::optional<core::NumericalRank> conditioning;
    double residual_norm = std::numeric_limits<double>::infinity();
    auto make_record = [&](double residual, std::size_t iterations,
                           core::Stop stop, bool classify) {
        core::Convergence result;
        result.residual = residual;
        if (conditioning) result.rank = std::make_shared<core::NumericalRank>(*conditioning);
        if (result.rank) if (auto& warning = result.rank->first_conditioning_warning) {
            warning->controls = {temperature, molar_density};
            warning->controls.insert(warning->controls.end(), mole_fractions.begin(), mole_fractions.end());
            if (core::AssociationObservations::current()) core::AssociationObservations::publish(*result.rank);
        }
        core::stopping_criterion(result, options, policy::association);
        result.floor = std::numeric_limits<double>::epsilon() * sites;
        result.amplification = last_pivot;
        result.contraction = history.contraction();
        result.iterations = iterations;
        result.evaluations = evaluations;
        result.stop = classify
            ? core::classify(residual, *result.criterion, result.floor,
                             result.contraction, true, result.amplification)
            : stop;
        return result;
    };
    const auto max_iterations = static_cast<std::size_t>(
        policy::ASSOCIATION_MAX_ITERATIONS.value);
    const double tolerance = options.value(policy::association);
    // Hybrid Newton / damped successive substitution (Kakalis et al.; cited by
    // Pereira2012 for this system). A damped Newton step is taken when some step
    // length inside the physical domain reduces the residual; otherwise the damped
    // substitution step is used. The converged point and the tolerance are
    // unchanged -- only the path to it.
    auto substitution_step = [&](const std::vector<double>& current) {
        auto step = bonded_densities(mole_fractions, parameters, current, delta);
        for (std::size_t a = 0; a < sites; ++a)
            step[a] = 0.5 * (current[a] + 1.0 / (1.0 + number_density * step[a]));
        return step;
    };
    auto residual_infinity_norm = [&](const std::vector<double>& point) {
        const auto values = association_residuals(
            molar_density, mole_fractions, parameters, point, delta
        );
        ++evaluations;
        double norm = 0.0;
        for (const double value : values) {
            const double magnitude = std::abs(value);
            if (!std::isfinite(magnitude)) return std::numeric_limits<double>::infinity();
            norm = std::max(norm, magnitude);
        }
        return norm;
    };
    auto in_domain = [&](const std::vector<double>& point) {
        for (const double value : point) {
            if (!std::isfinite(value)
                || value <= policy::SITE_FRACTION_FLOOR.value
                || value > policy::SITE_FRACTION_CEILING.value) {
                return false;
            }
        }
        return true;
    };

    for (std::size_t iteration = 0; iteration < max_iterations; ++iteration) {
        std::vector<double> next;
        last_pivot.reset();
        double damping = 0.0;
        std::optional<std::string> method;
        const bool detailed = scope.diagnostics_on();
        {
            const auto bonded = bonded_densities(mole_fractions, parameters, fractions, delta);
            std::vector<double> negative_residual(sites, 0.0);
            for (std::size_t a = 0; a < sites; ++a)
                negative_residual[a] = 1.0 - fractions[a] * (1.0 + number_density * bonded[a]);
            // Chapman1990 mass action: exact equivalent-partner pair, otherwise row quadratic.
            if (iteration == 0 && ((pair_start && start == nullptr) || number_density * *std::max_element(bonded.begin(), bonded.end())
                > 1.0 / policy::ASSOCIATION_JACOBIAN_CONDITION.value)) {
                next.resize(sites);
                for (std::size_t a = 0; a < sites; ++a) {
                    const double q = number_density * bonded[a];
                    double r = q;
                    bool found = false;
                    for (std::size_t b = 0; b < sites; ++b) {
                        if (mole_fractions[parameters.site_component[b]]
                            * parameters.site_multiplicity[b] * delta[a * sites + b] <= 0.0) continue;
                        const double partner = number_density * bonded[b];
                        if (found && r != partner) { r = q; break; }
                        r = partner;
                        found = true;
                    }
                    const double c = 1.0 + (q - r);
                    const double d = std::hypot(c, 2.0 * std::sqrt(r));
                    next[a] = c >= 0.0 ? 2.0 / (c + d) : (d - c) / (2.0 * r);
                }
                if (!in_domain(next)) next.clear();
                else if (detailed) scope.observe([&] { method.emplace("analytic association start"); });
            }
            double pivot_ratio = 0.0;
            const auto solved = next.empty() ? solve_linear(mass_action_jacobian(mole_fractions, parameters, fractions,
                number_density, delta, bonded), negative_residual, &pivot_ratio, &first, &conditioning, fractions, iteration)
                : std::nullopt;
            // singular or ill-conditioned: fall through to substitution
            if (next.empty() && !solved && detailed)
                scope.observe([&] { method.emplace("association Jacobian is ill-conditioned"); });
            if (solved) {
                const auto& step = *solved;
                // Damped Newton: cap the step inside (0, 1], then halve until the residual decreases.
                double alpha = 1.0;
                for (std::size_t a = 0; a < sites; ++a) {
                    if (step[a] < 0.0)
                        alpha = std::min(alpha, policy::ASSOCIATION_BOUNDARY_FRACTION.value * fractions[a] / -step[a]);
                    else if (step[a] > 0.0 && fractions[a] + step[a] > 1.0)
                        alpha = std::min(alpha, (1.0 - fractions[a]) / step[a]);
                }
                const double current = iteration == 0 ? residual_infinity_norm(fractions) : residual_norm;
                for (; next.empty() && alpha > std::numeric_limits<double>::epsilon(); alpha *= 0.5) {
                    std::vector<double> candidate(sites, 0.0);
                    for (std::size_t a = 0; a < sites; ++a) candidate[a] = fractions[a] + alpha * step[a];
                    if (in_domain(candidate) && residual_infinity_norm(candidate) < current) {
                        next = std::move(candidate);
                        last_pivot = pivot_ratio;
                        if (detailed) damping = alpha;
                    }
                }
                if (detailed && next.empty()) scope.observe([&] { method.emplace("Newton line search exhausted"); });
            }
        }
        if (next.empty()) {
            next = substitution_step(fractions);
        }
        const auto residuals = association_residuals(
            molar_density, mole_fractions, parameters, next, delta
        );
        const bool valid = in_domain(next);
        bool finite_residual = true;
        residual_norm = valid ? 0.0 : std::numeric_limits<double>::infinity();
        for (std::size_t site = 0; valid && site < sites; ++site) {
            const double magnitude = std::abs(residuals[site]);
            if (!std::isfinite(magnitude)) {
                finite_residual = false;
                residual_norm = std::numeric_limits<double>::infinity();
                break;
            }
            residual_norm = std::max(residual_norm, magnitude);
        }
        ++evaluations;
        if (detailed) scope.probe("association_iteration",
            {residual_norm, damping, last_pivot.value_or(std::numeric_limits<double>::quiet_NaN()),
             last_pivot ? 1.0 : 0.0}, static_cast<std::uint32_t>(iteration),
            last_pivot ? "Newton" : method ? method->c_str() : "");
        if (std::isfinite(residual_norm)) {
            history.push(residual_norm);
            if (!detailed) scope.probe("association_residual_norm", residual_norm,
                        static_cast<std::uint32_t>(iteration));
        }
        if (!valid) {
            const auto failure = make_record(
                residual_norm, iteration + 1, core::Stop::left_domain, false
            );
            scope.finish(failure);
            throw core::ConvergenceFailure(failure, first);
        }
        if (!finite_residual) {
            const auto failure = make_record(
                residual_norm, iteration + 1, core::Stop::nonfinite, false
            );
            scope.finish(failure);
            throw core::ConvergenceFailure(failure, first);
        }
        fractions = std::move(next);
        if (residual_norm <= tolerance) {
            const auto convergence = make_record(
                residual_norm, iteration + 1, core::Stop::converged, false
            );
            if (record != nullptr) {
                *record = convergence;
                if (record->amplification.has_value()) {
                    scope.probe("association_pivot_ratio", *record->amplification);
                }
                scope.finish(*record);
            }
            if (detailed && record == nullptr) scope.finish(convergence);
            scope.probe("site_fractions", fractions);
            return fractions;
        }
    }
    const auto failure = make_record(
        residual_norm, max_iterations, core::Stop::budget_exhausted, true
    );
    scope.finish(failure);
    throw core::ConvergenceFailure(failure, first);
}


}  // namespace epcsaft
