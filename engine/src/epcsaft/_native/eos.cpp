#include "eos.hpp"

#include "association.hpp"
#include "born.hpp"
#include "density.hpp"
#include "derivatives.hpp"
#include "electrolyte.hpp"
#include "helmholtz.hpp"
#include "phase_block.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace epcsaft {
namespace {

bool positive_finite(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value) && value > 0.0;
    });
}

bool finite(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
    });
}

}  // namespace

NativeModel::NativeModel(NativeParameters parameters)
    : parameters_(std::move(parameters)) {
    const std::size_t components = parameters_.segment_count.size();
    const std::size_t terms = parameters_.sigma_term_component.size();
    const std::size_t sites = parameters_.site_component.size();
    if (components == 0 || parameters_.sigma_constant.size() != components ||
        parameters_.dispersion_energy_over_k.size() != components ||
        parameters_.k_ij.size() != components * components ||
        parameters_.sigma_term_amplitude.size() != terms ||
        parameters_.sigma_term_exponent.size() != terms ||
        parameters_.site_multiplicity.size() != sites ||
        parameters_.association_active.size() != sites * sites ||
        parameters_.association_energy_over_k.size() != sites * sites ||
        parameters_.association_volume.size() != sites * sites ||
        parameters_.molar_mass_kg_per_mol.size() != components ||
        parameters_.charge_number.size() != components ||
        parameters_.relative_permittivity.size() != components ||
        parameters_.born_diameter_angstrom.size() != components ||
        parameters_.solvation_factor.size() != components) {
        throw std::invalid_argument("inconsistent native model dimensions");
    }
    if (!positive_finite(parameters_.segment_count) ||
        !positive_finite(parameters_.sigma_constant) ||
        !positive_finite(parameters_.dispersion_energy_over_k) ||
        !finite(parameters_.k_ij) ||
        !finite(parameters_.sigma_term_amplitude) ||
        !finite(parameters_.sigma_term_exponent) ||
        !finite(parameters_.association_energy_over_k) ||
        !finite(parameters_.association_volume) ||
        !finite(parameters_.molar_mass_kg_per_mol) ||
        !finite(parameters_.relative_permittivity) ||
        !std::isfinite(parameters_.dielectric_ion_suppression_coefficient) ||
        !finite(parameters_.born_diameter_angstrom) ||
        !finite(parameters_.solvation_factor) ||
        !std::isfinite(parameters_.ionic_region_relative_permittivity)) {
        throw std::invalid_argument("invalid native model parameters");
    }
    for (int component : parameters_.sigma_term_component) {
        if (component < 0 || static_cast<std::size_t>(component) >= components) {
            throw std::invalid_argument("invalid sigma correlation component");
        }
    }
    for (std::size_t site = 0; site < sites; ++site) {
        const int component = parameters_.site_component[site];
        if (component < 0 || static_cast<std::size_t>(component) >= components ||
            parameters_.site_multiplicity[site] <= 0) {
            throw std::invalid_argument("invalid association site");
        }
    }
    for (std::size_t pair = 0; pair < sites * sites; ++pair) {
        if (parameters_.association_active[pair] != 0 &&
            (parameters_.association_energy_over_k[pair] <= 0.0 ||
             parameters_.association_volume[pair] <= 0.0)) {
            throw std::invalid_argument("invalid association pair");
        }
    }
    const bool ionic = std::any_of(
        parameters_.charge_number.begin(),
        parameters_.charge_number.end(),
        [](int charge) { return charge != 0; }
    );
    if (ionic) {
        if (!(parameters_.dielectric_ion_suppression_coefficient > 0.0)
            || !(parameters_.ionic_region_relative_permittivity > 0.0)) {
            throw std::invalid_argument("invalid electrolyte model parameters");
        }
        for (std::size_t component = 0; component < components; ++component) {
            if (!(parameters_.solvation_factor[component] > 0.0)) {
                throw std::invalid_argument("invalid solvation factor");
            }
            if (parameters_.charge_number[component] == 0
                && (!(parameters_.molar_mass_kg_per_mol[component] > 0.0)
                    || !(parameters_.relative_permittivity[component] > 0.0))) {
                throw std::invalid_argument("invalid solvent dielectric parameters");
            }
            if (parameters_.charge_number[component] != 0
                && !(parameters_.born_diameter_angstrom[component] > 0.0)) {
                throw std::invalid_argument("invalid fitted Born diameter");
            }
        }
    }
}

bool NativeModel::supports_native_sdk_v1() const {
    return detail::supports_pure_phase_block(parameters_);
}

NativePurePhaseBlock NativeModel::evaluate_pure_phase(
    double temperature,
    double amount,
    double volume
) const {
    return detail::pure_phase_block(temperature, amount, volume, parameters_);
}

NativeParameterizedPurePhaseBlock NativeModel::evaluate_pure_phase_parameters(
    double temperature,
    double amount,
    double volume,
    double segment_count,
    double segment_diameter_angstrom,
    double dispersion_energy_over_k
) const {
    return detail::parameterized_pure_phase_block(
        temperature,
        amount,
        volume,
        segment_count,
        segment_diameter_angstrom,
        dispersion_energy_over_k,
        parameters_
    );
}

NativeResult NativeModel::evaluate(
    double temperature,
    double molar_density,
    const std::vector<double>& mole_fractions
) const {
    const std::size_t components = parameters_.segment_count.size();
    if (!std::isfinite(temperature) || temperature <= 0.0 ||
        !std::isfinite(molar_density) || molar_density <= 0.0 ||
        mole_fractions.size() != components || !finite(mole_fractions)) {
        throw std::invalid_argument("invalid native state");
    }
    const auto sigma = detail::segment_sigma(temperature, parameters_);
    if (!positive_finite(sigma)) {
        throw std::invalid_argument("invalid segment diameter correlation");
    }
    const auto diameter = detail::segment_diameter(temperature, parameters_, sigma);
    const auto site_fraction = detail::solve_site_fractions(
        temperature,
        molar_density,
        mole_fractions,
        parameters_,
        sigma,
        diameter
    );
    const auto electrolyte = detail::electrolyte_contribution(
        temperature,
        molar_density,
        temperature,
        molar_density,
        mole_fractions,
        mole_fractions,
        parameters_,
        sigma,
        sigma
    );
    const double born_ssm_ds = detail::born_ssm_ds_contribution(
        temperature,
        mole_fractions,
        parameters_,
        electrolyte.bulk_relative_permittivity
    );
    const auto values = detail::residual_helmholtz(
        temperature,
        molar_density,
        mole_fractions,
        parameters_,
        sigma,
        diameter,
        site_fraction,
        electrolyte.debye_huckel,
        born_ssm_ds
    );

    const double residual =
        values.hard_chain + values.dispersion + values.association
        + values.debye_huckel + values.born_ssm_ds;
    const auto derivatives = detail::thermodynamic_derivatives(
        temperature, molar_density, mole_fractions, parameters_, site_fraction
    );
    const double compressibility = 1.0 + molar_density * derivatives.density;
    double weighted_composition_derivative = 0.0;
    for (std::size_t component = 0; component < components; ++component) {
        weighted_composition_derivative += mole_fractions[component]
            * derivatives.composition[component];
    }
    std::vector<double> chemical_potential(components);
    std::optional<std::vector<double>> log_fugacity;
    if (std::isfinite(compressibility) && compressibility > 0.0) {
        log_fugacity.emplace(components);
    }
    for (std::size_t component = 0; component < components; ++component) {
        chemical_potential[component] = residual + compressibility - 1.0
            + derivatives.composition[component]
            - weighted_composition_derivative;
        if (log_fugacity.has_value()) {
            (*log_fugacity)[component] = chemical_potential[component]
                - std::log(compressibility);
        }
    }
    const double z_temperature = molar_density
        * derivatives.temperature_density;
    const double z_density = derivatives.density
        + molar_density * derivatives.density_density;
    const double pressure_slope_factor = compressibility
        + molar_density * z_density;
    const double cv_over_r = -2.0 * temperature * derivatives.temperature
        - temperature * temperature * derivatives.temperature_temperature;
    std::optional<double> residual_cp;
    if (std::isfinite(pressure_slope_factor)
        && pressure_slope_factor > 1.0e-10) {
        const double cp_over_r = cv_over_r
            + (compressibility + temperature * z_temperature)
                * (compressibility + temperature * z_temperature)
                / pressure_slope_factor
            - 1.0;
        if (std::isfinite(cp_over_r)) {
            residual_cp = detail::GAS_CONSTANT * cp_over_r;
        }
    }
    return {
        values.hard_chain,
        values.dispersion,
        values.association,
        values.debye_huckel,
        values.born_ssm_ds,
        residual,
        compressibility,
        compressibility * molar_density * detail::GAS_CONSTANT * temperature,
        std::move(chemical_potential),
        std::move(log_fugacity),
        detail::GAS_CONSTANT
            * (-residual - temperature * derivatives.temperature),
        -detail::GAS_CONSTANT * temperature * temperature
            * derivatives.temperature,
        detail::GAS_CONSTANT * temperature
            * (-temperature * derivatives.temperature + compressibility - 1.0),
        detail::GAS_CONSTANT * cv_over_r,
        residual_cp,
        molar_density,
        detail::GAS_CONSTANT * temperature * pressure_slope_factor,
        std::nullopt,
    };
}

NativeResult NativeModel::evaluate_pressure(
    double temperature,
    double pressure,
    const std::vector<double>& mole_fractions,
    const std::string& phase
) const {
    const std::size_t components = parameters_.segment_count.size();
    if (!std::isfinite(temperature) || temperature <= 0.0
        || !std::isfinite(pressure) || pressure <= 0.0
        || mole_fractions.size() != components || !finite(mole_fractions)) {
        throw std::invalid_argument("invalid native pressure state");
    }
    const auto sigma = detail::segment_sigma(temperature, parameters_);
    const auto diameter = detail::segment_diameter(temperature, parameters_, sigma);
    double packing_coefficient = 0.0;
    for (std::size_t component = 0; component < components; ++component) {
        packing_coefficient += mole_fractions[component]
            * parameters_.segment_count[component]
            * diameter[component] * diameter[component] * diameter[component];
    }
    packing_coefficient *= detail::PI / 6.0 * (detail::AVOGADRO / 1.0e30);
    if (!(packing_coefficient > 0.0) || !std::isfinite(packing_coefficient)) {
        throw std::invalid_argument("invalid pressure-state packing coefficient");
    }
    const double close_packing = detail::PI / (3.0 * std::sqrt(2.0));
    const auto root = detail::solve_density(
        pressure,
        pressure / (detail::GAS_CONSTANT * temperature),
        close_packing / packing_coefficient,
        phase,
        [this, temperature, &mole_fractions](double density) {
            return evaluate(temperature, density, mole_fractions);
        }
    );
    auto result = evaluate(temperature, root.molar_density, mole_fractions);
    result.density_diagnostics = NativeDensityDiagnostics{
        root.bracket_lower,
        root.bracket_upper,
        root.pressure_residual,
        root.pressure_density_derivative,
        root.iterations,
        root.stable,
        root.branch,
    };
    return result;
}

}  // namespace epcsaft
