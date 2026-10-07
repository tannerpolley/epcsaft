#pragma once

#include "core/convergence.hpp"
#include "core/numerical_policy.hpp"
#include "core/parameters.hpp"
#include "core/physical_constants.hpp"
#include "core/scalar.hpp"
#include "core/trace.hpp"
#include "contributions/helmholtz.hpp"

#include <cppad/cppad.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace epcsaft {

template <typename Scalar, typename ParameterSet>
std::vector<Scalar> association_delta_values(
    const Scalar& temperature,
    const Scalar& molar_density,
    const std::vector<Scalar>& mole_fractions,
    const ParameterSet& parameters,
    const std::vector<Scalar>& sigma,
    const std::vector<Scalar>& diameter,
    const std::vector<Scalar>& association_energy_over_k,
    const std::vector<Scalar>& association_volume
) {
    const std::size_t components = parameters.segment_count.size();
    const std::size_t sites = parameters.site_component.size();
    const auto zeta = packing_moments(molar_density, mole_fractions, parameters, diameter);
    const Scalar one_minus_eta = Scalar(1.0) - zeta[3];

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
            const Scalar contact = hard_sphere_contact(
                diameter[i] * diameter[j] / (diameter[i] + diameter[j]), zeta[2], one_minus_eta);
            // Wolbach--Sandler pairs (active == 2): the Gross2002 Eq. 3 size
            // factor times sigma_ij^3 is the geometric-mean diameter cubed.
            const Scalar size_ij = (parameters.association_active[pair] == 2
                    ? scalar_sqrt(sigma[i] * sigma[j]) : Scalar(0.5) * (sigma[i] + sigma[j]))
                * (Scalar(1.0) - parameters.l_ij[i * components + j]);
            // EqID: association_strength
            delta[pair] = contact
                * (scalar_exp(
                       association_energy_over_k[pair] / temperature
                   ) - Scalar(1.0))
                * size_ij * size_ij * size_ij
                * association_volume[pair];
        }
    }
    return delta;
}

template <typename Scalar>
// EqExclude: scalar utility dispatch
std::vector<Scalar> association_delta(
    const Scalar& temperature,
    const Scalar& molar_density,
    const std::vector<Scalar>& mole_fractions,
    const Parameters& parameters,
    const std::vector<Scalar>& sigma,
    const std::vector<Scalar>& diameter
) {
    return association_delta_values(
        temperature,
        molar_density,
        mole_fractions,
        parameters,
        sigma,
        diameter,
        std::vector<Scalar>(
            parameters.association_energy_over_k.begin(),
            parameters.association_energy_over_k.end()
        ),
        std::vector<Scalar>(
            parameters.association_volume.begin(),
            parameters.association_volume.end()
        )
    );
}

template <typename Scalar, typename ParameterSet>
std::vector<Scalar> association_residuals(
    const Scalar& molar_density,
    const std::vector<Scalar>& mole_fractions,
    const ParameterSet& parameters,
    const std::vector<Scalar>& site_fractions,
    const std::vector<Scalar>& delta
) {
    const std::size_t sites = parameters.site_component.size();
    const Scalar number_density =
        molar_density * NUMBER_DENSITY_PER_ANGSTROM3;
    std::vector<Scalar> residuals(sites, Scalar(0.0));
    for (std::size_t a = 0; a < sites; ++a) {
        Scalar bonded_density = Scalar(0.0);
        for (std::size_t b = 0; b < sites; ++b) {
            const auto component =
                static_cast<std::size_t>(parameters.site_component[b]);
            bonded_density += mole_fractions[component]
                * parameters.site_multiplicity[b]
                * site_fractions[b]
                * delta[a * sites + b];
        }
        // EqID: association_mass_action
        residuals[a] = site_fractions[a]
            * (Scalar(1.0) + number_density * bonded_density)
            - Scalar(1.0);
    }
    return residuals;
}

// EqExclude: numerical iteration mechanics
// Dense and small: the site count is single digits, so the factorisation never
// shows up next to the tape sweeps. Full pivoting because a nearly singular
// Jacobian is the signal to fall back to substitution, and the diagonal ratio
// reports it.
inline std::vector<double> solve_linear(
    const std::vector<double>& matrix,
    const std::vector<double>& right_hand_side,
    double* pivot_ratio = nullptr
) {
    using RowMajor =
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    const auto size = static_cast<Eigen::Index>(right_hand_side.size());
    const Eigen::FullPivLU<Eigen::MatrixXd> factorization(
        Eigen::Map<const RowMajor>(matrix.data(), size, size)
    );
    const auto pivots = factorization.matrixLU().diagonal().cwiseAbs();
    const double largest = pivots.maxCoeff();
    const double smallest = pivots.minCoeff();
    if (!std::isfinite(largest) || !(largest > 0.0)
        || smallest <= policy::ASSOCIATION_JACOBIAN_CONDITION.value * largest) {
        throw std::runtime_error("association Jacobian is ill-conditioned");
    }
    if (pivot_ratio != nullptr) *pivot_ratio = largest / smallest;
    const Eigen::VectorXd solution = factorization.solve(
        Eigen::Map<const Eigen::VectorXd>(right_hand_side.data(), size)
    );
    return {solution.data(), solution.data() + size};
}

// EqExclude: numerical iteration mechanics
inline std::vector<double> solve_site_fractions(
    double temperature,
    double molar_density,
    const std::vector<double>& mole_fractions,
    const Parameters& parameters,
    const std::vector<double>& sigma,
    const std::vector<double>& diameter,
    core::Convergence* record = nullptr
) {
    const std::size_t sites = parameters.site_component.size();
    if (sites == 0) {
        if (record != nullptr) {
            *record = {};
            record->criterion = &policy::ASSOCIATION_RESIDUAL;
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
    const auto delta = association_delta_values(
        temperature,
        molar_density,
        mole_fractions,
        parameters,
        sigma,
        diameter,
        parameters.association_energy_over_k,
        parameters.association_volume
    );
    const double number_density =
        molar_density * NUMBER_DENSITY_PER_ANGSTROM3;
    std::vector<double> fractions(sites, 1.0);
    core::ResidualHistory history;
    std::size_t evaluations = 0;
    std::optional<double> last_pivot;
    double residual_norm = std::numeric_limits<double>::infinity();
    auto make_record = [&](double residual, std::size_t iterations,
                           core::Stop stop, bool classify) {
        core::Convergence result;
        result.residual = residual;
        result.criterion = &policy::ASSOCIATION_RESIDUAL;
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
    const double tolerance = policy::ASSOCIATION_RESIDUAL.value;
    // Hybrid Newton / damped successive substitution (Kakalis et al.; cited by
    // Pereira2012 for this system). A damped Newton step is taken when some step
    // length inside the physical domain reduces the residual; otherwise the damped
    // substitution step is used. The converged point and the tolerance are
    // unchanged -- only the path to it.
    auto substitution_step = [&](const std::vector<double>& current) {
        std::vector<double> step(sites, 1.0);
        for (std::size_t a = 0; a < sites; ++a) {
            double bonded_density = 0.0;
            for (std::size_t b = 0; b < sites; ++b) {
                const auto component =
                    static_cast<std::size_t>(parameters.site_component[b]);
                bonded_density += mole_fractions[component]
                    * parameters.site_multiplicity[b]
                    * current[b]
                    * delta[a * sites + b];
            }
            const double fixed_point = 1.0 / (1.0 + number_density * bonded_density);
            step[a] = 0.5 * (current[a] + fixed_point);
        }
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
        // Newton on F_a = X_a (1 + rho_N B_a) - 1, with
        // dF_a/dX_b = X_a rho_N x_c(b) n_b delta_ab + kron_ab (1 + rho_N B_a).
        {
            std::vector<double> bonded(sites, 0.0);
            for (std::size_t a = 0; a < sites; ++a) {
                for (std::size_t b = 0; b < sites; ++b) {
                    const auto component =
                        static_cast<std::size_t>(parameters.site_component[b]);
                    bonded[a] += mole_fractions[component]
                        * parameters.site_multiplicity[b]
                        * fractions[b]
                        * delta[a * sites + b];
                }
            }
            std::vector<double> jacobian(sites * sites, 0.0);
            std::vector<double> negative_residual(sites, 0.0);
            for (std::size_t a = 0; a < sites; ++a) {
                for (std::size_t b = 0; b < sites; ++b) {
                    const auto component =
                        static_cast<std::size_t>(parameters.site_component[b]);
                    jacobian[a * sites + b] = fractions[a] * number_density
                        * mole_fractions[component]
                        * parameters.site_multiplicity[b]
                        * delta[a * sites + b];
                }
                jacobian[a * sites + a] += 1.0 + number_density * bonded[a];
                negative_residual[a] =
                    1.0 - fractions[a] * (1.0 + number_density * bonded[a]);
            }
            // Chapman1990 mass action: exact equivalent-partner pair, otherwise row quadratic.
            if (iteration == 0 && number_density * *std::max_element(bonded.begin(), bonded.end())
                > 1.0 / policy::ASSOCIATION_JACOBIAN_CONDITION.value) {
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
            if (next.empty()) try {
                double pivot_ratio = 0.0;
                const auto step = solve_linear(jacobian, negative_residual, &pivot_ratio);
                // Damped Newton: cap the step inside (0, 1], then halve until the residual decreases.
                double alpha = 1.0;
                for (std::size_t a = 0; a < sites; ++a) {
                    if (step[a] < 0.0)
                        alpha = std::min(alpha, policy::ASSOCIATION_BOUNDARY_FRACTION.value * fractions[a] / -step[a]);
                    else if (step[a] > 0.0)
                        alpha = std::min(alpha, (1.0 - fractions[a]) / step[a]);
                }
                const double current = residual_infinity_norm(fractions);
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
            } catch (const std::exception& error) {
                if (detailed) scope.observe([&] { method.emplace(error.what()); });
                // singular or ill-conditioned: fall through to substitution
            }
        }
        if (next.empty()) {
            next = substitution_step(fractions);
        }
        const auto residuals = association_residuals(
            molar_density, mole_fractions, parameters, next, delta
        );
        residual_norm = 0.0;
        bool valid = true;
        bool finite_residual = true;
        for (std::size_t site = 0; site < sites; ++site) {
            if (!std::isfinite(next[site])
                || next[site] <= policy::SITE_FRACTION_FLOOR.value
                || next[site] > policy::SITE_FRACTION_CEILING.value) {
                valid = false;
                residual_norm = std::numeric_limits<double>::infinity();
                break;
            }
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
            throw core::ConvergenceFailure(failure);
        }
        if (!finite_residual) {
            const auto failure = make_record(
                residual_norm, iteration + 1, core::Stop::nonfinite, false
            );
            scope.finish(failure);
            throw core::ConvergenceFailure(failure);
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
    throw core::ConvergenceFailure(failure);
}


}  // namespace epcsaft
