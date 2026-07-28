#pragma once

#include "association.hpp"
#include "born.hpp"
#include "electrolyte.hpp"
#include "helmholtz.hpp"
#include "model.hpp"

#include <cppad/cppad.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <exception>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

namespace epcsaft::detail {

namespace cppad_parallel {

inline std::recursive_mutex execution_mutex;
inline std::atomic<bool> active{false};
inline thread_local std::size_t thread_index = 0;
inline thread_local bool worker = false;

// EqExclude: transport and buffer mechanics
inline bool in_parallel() {
    return active.load(std::memory_order_acquire);
}

// EqExclude: transport and buffer mechanics
inline std::size_t thread_number() {
    return thread_index;
}

}  // namespace cppad_parallel

class CppAdSerialGuard final {
public:
    // EqExclude: transport and buffer mechanics
    CppAdSerialGuard()
        : lock_(cppad_parallel::execution_mutex, std::defer_lock) {
        if (!cppad_parallel::worker) {
            lock_.lock();
        }
    }

private:
    std::unique_lock<std::recursive_mutex> lock_;
};

template <typename Function>
// EqExclude: transport and buffer mechanics
void cppad_parallel_for(std::size_t count, Function&& function) {
    if (count == 0) {
        return;
    }
    const unsigned int hardware_threads = std::thread::hardware_concurrency();
    const std::size_t physical_thread_estimate = std::max<std::size_t>(
        1, hardware_threads == 0 ? 1 : hardware_threads / 2
    );
    const std::size_t thread_count = std::min(
        {count,
         physical_thread_estimate,
         static_cast<std::size_t>(CPPAD_MAX_NUM_THREADS)}
    );
    if (thread_count == 1) {
        CppAdSerialGuard guard;
        for (std::size_t index = 0; index < count; ++index) {
            function(index);
        }
        return;
    }

    std::unique_lock<std::recursive_mutex> execution_lock(
        cppad_parallel::execution_mutex
    );
    CppAD::thread_alloc::parallel_setup(
        thread_count,
        cppad_parallel::in_parallel,
        cppad_parallel::thread_number
    );
    CppAD::thread_alloc::hold_memory(true);
    CppAD::parallel_ad<double>();

    std::atomic<std::size_t> next_index{0};
    std::mutex exception_mutex;
    std::exception_ptr first_exception;
    cppad_parallel::active.store(true, std::memory_order_release);
    auto run = [&](std::size_t worker_index) {
        cppad_parallel::thread_index = worker_index;
        cppad_parallel::worker = true;
        try {
            while (true) {
                const std::size_t index = next_index.fetch_add(
                    1, std::memory_order_relaxed
                );
                if (index >= count) {
                    break;
                }
                function(index);
            }
        } catch (...) {
            std::lock_guard<std::mutex> exception_lock(exception_mutex);
            if (first_exception == nullptr) {
                first_exception = std::current_exception();
            }
            next_index.store(count, std::memory_order_relaxed);
        }
        cppad_parallel::worker = false;
        cppad_parallel::thread_index = 0;
    };

    auto reset_cppad = [] {
        cppad_parallel::active.store(false, std::memory_order_release);
        CppAD::thread_alloc::parallel_setup(1, nullptr, nullptr);
        CppAD::thread_alloc::hold_memory(false);
        CppAD::parallel_ad<double>();
    };
    std::vector<std::thread> workers;
    try {
        workers.reserve(thread_count - 1);
        for (std::size_t index = 1; index < thread_count; ++index) {
            workers.emplace_back(run, index);
        }
        run(0);
        for (auto& thread : workers) {
            thread.join();
        }
    } catch (...) {
        next_index.store(count, std::memory_order_relaxed);
        for (auto& thread : workers) {
            if (thread.joinable()) {
                thread.join();
            }
        }
        reset_cppad();
        throw;
    }
    reset_cppad();
    if (first_exception != nullptr) {
        std::rethrow_exception(first_exception);
    }
}

struct ThermodynamicDerivatives final {
    double temperature;
    double density;
    std::vector<double> composition;
    double temperature_temperature;
    double temperature_density;
    double density_density;
    double residual_helmholtz;
    std::optional<double> active_born;
    std::optional<double> active_solvation_factor;
    std::optional<double> active_ionic_permittivity;
    std::vector<double> active_aqueous_k_ij;
    std::optional<double> active_association_parameter;
    std::vector<double> hessian;
};

// EqExclude: CppAD-generated derivative machinery
inline double directional_second(
    const std::vector<double>& hessian,
    const std::vector<double>& left,
    const std::vector<double>& right
) {
    const std::size_t size = left.size();
    double value = 0.0;
    for (std::size_t row = 0; row < size; ++row) {
        for (std::size_t column = 0; column < size; ++column) {
            value += left[row] * hessian[row * size + column] * right[column];
        }
    }
    return value;
}

// EqExclude: CppAD-generated derivative machinery
inline ThermodynamicDerivatives thermodynamic_derivatives_impl(
    double temperature,
    double molar_density,
    const std::vector<double>& mole_fractions,
    const NativeParameters& parameters,
    const std::vector<double>& site_fractions,
    std::optional<std::size_t> active_born_component,
    double active_born_diameter_angstrom,
    const std::optional<double>& active_solvation_factor,
    const std::optional<double>& active_ionic_permittivity,
    const std::optional<std::array<double, 3>>& active_aqueous_k_ij,
    const std::optional<double>& active_association_parameter,
    bool active_association_energy
) {
    CppAdSerialGuard cppad_guard;
    using ADScalar = CppAD::AD<double>;
    const std::size_t components = mole_fractions.size();
    const std::size_t sites = site_fractions.size();
    const std::size_t born_coordinates =
        static_cast<std::size_t>(active_born_component.has_value());
    const std::size_t solvation_coordinates =
        static_cast<std::size_t>(active_solvation_factor.has_value());
    const std::size_t ionic_permittivity_coordinates =
        static_cast<std::size_t>(active_ionic_permittivity.has_value());
    const std::size_t kij_coordinates = active_aqueous_k_ij.has_value() ? 3 : 0;
    const std::size_t association_coordinates =
        static_cast<std::size_t>(active_association_parameter.has_value());
    const std::size_t coordinates = 2 + components + born_coordinates
        + solvation_coordinates + ionic_permittivity_coordinates
        + kij_coordinates + association_coordinates;
    const std::size_t born_offset = 2 + components;
    const std::size_t solvation_offset = born_offset + born_coordinates;
    const std::size_t ionic_permittivity_offset =
        solvation_offset + solvation_coordinates;
    const std::size_t kij_offset =
        ionic_permittivity_offset + ionic_permittivity_coordinates;
    const std::size_t association_offset = kij_offset + kij_coordinates;
    const std::size_t width = coordinates + sites;

    std::vector<ADScalar> independent(width);
    independent[0] = temperature;
    independent[1] = molar_density;
    for (std::size_t component = 0; component < components; ++component) {
        independent[2 + component] = mole_fractions[component];
    }
    if (active_born_component.has_value()) {
        independent[born_offset] = active_born_diameter_angstrom;
    }
    if (active_solvation_factor.has_value()) {
        independent[solvation_offset] = *active_solvation_factor;
    }
    if (active_ionic_permittivity.has_value()) {
        independent[ionic_permittivity_offset] = *active_ionic_permittivity;
    }
    if (active_aqueous_k_ij.has_value()) {
        for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
            independent[kij_offset + coordinate] = (*active_aqueous_k_ij)[coordinate];
        }
    }
    if (active_association_parameter.has_value()) {
        independent[association_offset] = *active_association_parameter;
    }
    for (std::size_t site = 0; site < sites; ++site) {
        independent[coordinates + site] = site_fractions[site];
    }
    CppAD::Independent(independent);

    const std::vector<ADScalar> composition(
        independent.begin() + 2, independent.begin() + 2 + components
    );
    const std::vector<ADScalar> fractions(
        independent.begin() + coordinates, independent.end()
    );
    const auto sigma = segment_sigma(independent[0], parameters);
    const auto diameter = segment_diameter(independent[0], parameters, sigma);
    const auto reference_sigma = segment_sigma(temperature, parameters);
    const auto reference_diameter = segment_diameter(
        temperature, parameters, reference_sigma
    );
    const auto electrolyte = electrolyte_contribution(
        independent[0],
        independent[1],
        temperature,
        molar_density,
        composition,
        mole_fractions,
        parameters,
        diameter,
        reference_diameter
    );
    const auto born = active_born_component.has_value()
        ? born_ssm_ds_contribution_active(
            independent[0],
            composition,
            parameters,
            electrolyte.bulk_relative_permittivity,
            *active_born_component,
            independent[born_offset]
        )
        : active_solvation_factor.has_value()
            ? born_ssm_ds_contribution_active_solvation_factor(
                independent[0],
                composition,
                parameters,
                electrolyte.bulk_relative_permittivity,
                0,
                independent[solvation_offset]
            )
            : active_ionic_permittivity.has_value()
                ? born_ssm_ds_contribution_active_ionic_permittivity(
                    independent[0],
                    composition,
                    parameters,
                    electrolyte.bulk_relative_permittivity,
                    independent[ionic_permittivity_offset]
                )
            : born_contribution(
            independent[0],
            composition,
            parameters,
            electrolyte.bulk_relative_permittivity,
            sigma
        );
    std::vector<ADScalar> active_k_values;
    if (active_aqueous_k_ij.has_value()) {
        active_k_values.assign(
            independent.begin() + static_cast<std::ptrdiff_t>(kij_offset),
            independent.begin() + static_cast<std::ptrdiff_t>(kij_offset + 3)
        );
    }
    const auto contributions = residual_helmholtz(
        independent[0],
        independent[1],
        composition,
        parameters,
        sigma,
        diameter,
        fractions,
        electrolyte.debye_huckel,
        born,
        active_k_values
    );
    const auto delta = active_association_parameter.has_value()
        ? association_delta_active_pure_parameter(
            independent[0],
            independent[1],
            composition,
            parameters,
            sigma,
            diameter,
            independent[association_offset],
            active_association_energy
        )
        : association_delta(
            independent[0],
            independent[1],
            composition,
            parameters,
            sigma,
            diameter
        );
    const auto mass_action = association_residuals(
        independent[1], composition, parameters, fractions, delta
    );
    std::vector<ADScalar> dependent = {
        contributions.hard_chain + contributions.dispersion
            + contributions.association + contributions.debye_huckel
            + contributions.born_ssm_ds,
    };
    dependent.insert(dependent.end(), mass_action.begin(), mass_action.end());
    CppAD::ADFun<double> function(independent, dependent);

    std::vector<double> point = {temperature, molar_density};
    point.insert(point.end(), mole_fractions.begin(), mole_fractions.end());
    if (active_born_component.has_value()) {
        point.push_back(active_born_diameter_angstrom);
    }
    if (active_solvation_factor.has_value()) {
        point.push_back(*active_solvation_factor);
    }
    if (active_ionic_permittivity.has_value()) {
        point.push_back(*active_ionic_permittivity);
    }
    if (active_aqueous_k_ij.has_value()) {
        point.insert(
            point.end(), active_aqueous_k_ij->begin(), active_aqueous_k_ij->end()
        );
    }
    if (active_association_parameter.has_value()) {
        point.push_back(*active_association_parameter);
    }
    point.insert(point.end(), site_fractions.begin(), site_fractions.end());
    const auto jacobian = function.Jacobian(point);

    std::vector<double> association_jacobian(sites * sites);
    for (std::size_t row = 0; row < sites; ++row) {
        for (std::size_t column = 0; column < sites; ++column) {
            association_jacobian[row * sites + column] =
                jacobian[(1 + row) * width + coordinates + column];
        }
    }
    std::vector<std::vector<double>> site_sensitivity(
        coordinates, std::vector<double>(sites, 0.0)
    );
    for (std::size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
        if (sites == 0) {
            break;
        }
        std::vector<double> right_hand_side(sites);
        for (std::size_t row = 0; row < sites; ++row) {
            right_hand_side[row] = -jacobian[(1 + row) * width + coordinate];
        }
        site_sensitivity[coordinate] = solve_linear(
            association_jacobian, std::move(right_hand_side)
        );
    }

    std::vector<double> gradient(coordinates);
    for (std::size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
        gradient[coordinate] = jacobian[coordinate];
        for (std::size_t site = 0; site < sites; ++site) {
            gradient[coordinate] += jacobian[coordinates + site]
                * site_sensitivity[coordinate][site];
        }
    }

    std::vector<std::vector<double>> directions(
        coordinates, std::vector<double>(width, 0.0)
    );
    for (std::size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
        directions[coordinate][coordinate] = 1.0;
        for (std::size_t site = 0; site < sites; ++site) {
            directions[coordinate][coordinates + site] =
                site_sensitivity[coordinate][site];
        }
    }
    const auto helmholtz_hessian = function.Hessian(point, 0);

    std::vector<std::vector<double>> association_hessians(sites);
    for (std::size_t site = 0; site < sites; ++site) {
        association_hessians[site] = function.Hessian(point, 1 + site);
    }
    auto reduced_second = [&](std::size_t left, std::size_t right) {
        std::vector<double> site_second(sites, 0.0);
        if (sites != 0) {
            std::vector<double> right_hand_side(sites);
            for (std::size_t site = 0; site < sites; ++site) {
                right_hand_side[site] = -directional_second(
                    association_hessians[site], directions[left], directions[right]
                );
            }
            site_second = solve_linear(
                association_jacobian, std::move(right_hand_side)
            );
        }
        double value = directional_second(
            helmholtz_hessian, directions[left], directions[right]
        );
        for (std::size_t site = 0; site < sites; ++site) {
            value += jacobian[coordinates + site] * site_second[site];
        }
        return value;
    };

    std::vector<double> reduced_hessian(coordinates * coordinates);
    for (std::size_t row = 0; row < coordinates; ++row) {
        for (std::size_t column = 0; column < coordinates; ++column) {
            reduced_hessian[row * coordinates + column] =
                reduced_second(row, column);
        }
    }
    const double residual_helmholtz = function.Forward(0, point)[0];

    return {
        gradient[0],
        gradient[1],
        std::vector<double>(
            gradient.begin() + 2,
            gradient.begin() + 2 + static_cast<std::ptrdiff_t>(components)
        ),
        reduced_hessian[0],
        reduced_hessian[1],
        reduced_hessian[coordinates + 1],
        residual_helmholtz,
        active_born_component.has_value()
            ? std::optional<double>(gradient[born_offset])
            : std::nullopt,
        active_solvation_factor.has_value()
            ? std::optional<double>(gradient[solvation_offset])
            : std::nullopt,
        active_ionic_permittivity.has_value()
            ? std::optional<double>(gradient[ionic_permittivity_offset])
            : std::nullopt,
        active_aqueous_k_ij.has_value()
            ? std::vector<double>(
                gradient.begin() + static_cast<std::ptrdiff_t>(kij_offset),
                gradient.begin() + static_cast<std::ptrdiff_t>(kij_offset + 3)
            )
            : std::vector<double>{},
        active_association_parameter.has_value()
            ? std::optional<double>(gradient[association_offset])
            : std::nullopt,
        std::move(reduced_hessian),
    };
}

// EqExclude: CppAD-generated derivative machinery
inline ThermodynamicDerivatives thermodynamic_derivatives(
    double temperature,
    double molar_density,
    const std::vector<double>& mole_fractions,
    const NativeParameters& parameters,
    const std::vector<double>& site_fractions
) {
    return thermodynamic_derivatives_impl(
        temperature,
        molar_density,
        mole_fractions,
        parameters,
        site_fractions,
        std::nullopt,
        0.0,
        std::nullopt,
        std::nullopt,
        std::nullopt,
        std::nullopt,
        false
    );
}

// EqExclude: CppAD-generated derivative machinery
inline ThermodynamicDerivatives thermodynamic_derivatives_active_born(
    double temperature,
    double molar_density,
    const std::vector<double>& mole_fractions,
    const NativeParameters& parameters,
    const std::vector<double>& site_fractions,
    std::size_t active_component,
    double active_born_diameter_angstrom
) {
    return thermodynamic_derivatives_impl(
        temperature,
        molar_density,
        mole_fractions,
        parameters,
        site_fractions,
        active_component,
        active_born_diameter_angstrom,
        std::nullopt,
        std::nullopt,
        std::nullopt,
        std::nullopt,
        false
    );
}

// EqExclude: CppAD-generated derivative machinery
inline ThermodynamicDerivatives thermodynamic_derivatives_active_solvation_factor(
    double temperature,
    double molar_density,
    const std::vector<double>& mole_fractions,
    const NativeParameters& parameters,
    const std::vector<double>& site_fractions,
    double active_solvation_factor
) {
    return thermodynamic_derivatives_impl(
        temperature,
        molar_density,
        mole_fractions,
        parameters,
        site_fractions,
        std::nullopt,
        0.0,
        active_solvation_factor,
        std::nullopt,
        std::nullopt,
        std::nullopt,
        false
    );
}

// EqExclude: CppAD-generated derivative machinery
inline ThermodynamicDerivatives
thermodynamic_derivatives_active_ionic_region_permittivity(
    double temperature,
    double molar_density,
    const std::vector<double>& mole_fractions,
    const NativeParameters& parameters,
    const std::vector<double>& site_fractions,
    double active_ionic_region_relative_permittivity
) {
    return thermodynamic_derivatives_impl(
        temperature,
        molar_density,
        mole_fractions,
        parameters,
        site_fractions,
        std::nullopt,
        0.0,
        std::nullopt,
        active_ionic_region_relative_permittivity,
        std::nullopt,
        std::nullopt,
        false
    );
}

// EqExclude: CppAD-generated derivative machinery
inline ThermodynamicDerivatives thermodynamic_derivatives_active_aqueous_kij(
    double temperature,
    double molar_density,
    const std::vector<double>& mole_fractions,
    const NativeParameters& parameters,
    const std::vector<double>& site_fractions,
    const std::array<double, 3>& active_aqueous_k_ij
) {
    return thermodynamic_derivatives_impl(
        temperature,
        molar_density,
        mole_fractions,
        parameters,
        site_fractions,
        std::nullopt,
        0.0,
        std::nullopt,
        std::nullopt,
        active_aqueous_k_ij,
        std::nullopt,
        false
    );
}

// EqExclude: CppAD-generated derivative machinery
inline ThermodynamicDerivatives thermodynamic_derivatives_active_pure_association_parameter(
    double temperature,
    double molar_density,
    const std::vector<double>& mole_fractions,
    const NativeParameters& parameters,
    const std::vector<double>& site_fractions,
    double active_parameter,
    bool active_energy
) {
    return thermodynamic_derivatives_impl(
        temperature,
        molar_density,
        mole_fractions,
        parameters,
        site_fractions,
        std::nullopt,
        0.0,
        std::nullopt,
        std::nullopt,
        std::nullopt,
        active_parameter,
        active_energy
    );
}

}  // namespace epcsaft::detail
