#include "eos.hpp"

#include "association.hpp"
#include "born.hpp"
#include "density.hpp"
#include "derivatives.hpp"
#include "electrolyte_reference.hpp"
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

detail::ReferenceObservation active_born_observation(
    const NativeModel& model,
    const NativeParameters& parameters,
    double temperature,
    double pressure,
    const std::vector<double>& composition,
    std::size_t active_component
) {
    const auto state = model.evaluate_pressure(
        temperature, pressure, composition, "liquid"
    );
    if (!state.log_fugacity_coefficient.has_value()) {
        throw std::runtime_error("reference fugacity is unavailable");
    }
    const auto sigma = detail::segment_sigma(temperature, parameters);
    const auto diameter = detail::segment_diameter(temperature, parameters, sigma);
    const auto site_fractions = detail::solve_site_fractions(
        temperature,
        state.molar_density,
        composition,
        parameters,
        sigma,
        diameter
    );
    const auto derivatives = detail::thermodynamic_derivatives_active_born(
        temperature,
        state.molar_density,
        composition,
        parameters,
        site_fractions,
        active_component,
        parameters.born_diameter_angstrom[active_component]
    );
    if (!derivatives.active_born.has_value()) {
        throw std::runtime_error("active Born derivative is unavailable");
    }
    const std::size_t components = composition.size();
    const std::size_t coordinates = 3 + components;
    const std::size_t active_coordinate = 2 + components;
    const auto hessian = [&derivatives, coordinates](
        std::size_t row, std::size_t column
    ) {
        return derivatives.hessian[row * coordinates + column];
    };
    const double compressibility = state.compressibility_factor;
    const double pressure_diameter = state.molar_density * state.molar_density
        * detail::GAS_CONSTANT * temperature * hessian(1, active_coordinate);
    const double density_diameter = -pressure_diameter
        / state.pressure_density_derivative;
    std::array<double, 2> total_derivative{};
    for (std::size_t ion = 0; ion < 2; ++ion) {
        const std::size_t component = 1 + ion;
        double weighted_composition_diameter = 0.0;
        double weighted_composition_density = 0.0;
        for (std::size_t mixture_component = 0;
             mixture_component < components;
             ++mixture_component) {
            weighted_composition_diameter += composition[mixture_component]
                * hessian(2 + mixture_component, active_coordinate);
            weighted_composition_density += composition[mixture_component]
                * hessian(2 + mixture_component, 1);
        }
        const double chemical_potential_diameter = *derivatives.active_born
            + state.molar_density * hessian(1, active_coordinate)
            + hessian(2 + component, active_coordinate)
            - weighted_composition_diameter;
        const double log_fugacity_diameter = chemical_potential_diameter
            - state.molar_density * hessian(1, active_coordinate)
                / compressibility;
        const double compressibility_density = derivatives.density
            + state.molar_density * derivatives.density_density;
        const double chemical_potential_density = 2.0 * derivatives.density
            + state.molar_density * derivatives.density_density
            + hessian(2 + component, 1)
            - weighted_composition_density;
        const double log_fugacity_density = chemical_potential_density
            - compressibility_density / compressibility;
        total_derivative[ion] = detail::GAS_CONSTANT * temperature
            * (log_fugacity_diameter
               + log_fugacity_density * density_diameter);
    }
    return {
        {
            (*state.log_fugacity_coefficient)[1],
            (*state.log_fugacity_coefficient)[2],
        },
        total_derivative,
    };
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
        parameters_.l_ij.size() != components * components ||
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
        !finite(parameters_.l_ij) ||
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
        if (parameters_.electrolyte_model != ELECTROLYTE_MODEL_FIGIEL_SSM_DS
            && parameters_.electrolyte_model
                != ELECTROLYTE_MODEL_ASCANI_ORIGINAL_BORN
            && parameters_.electrolyte_model
                != ELECTROLYTE_MODEL_HELD_2008_DEBYE_HUCKEL) {
            throw std::invalid_argument("invalid electrolyte model selection");
        }
        for (std::size_t component = 0; component < components; ++component) {
            if (parameters_.charge_number[component] == 0
                && (!(parameters_.molar_mass_kg_per_mol[component] > 0.0)
                    || !(parameters_.relative_permittivity[component] > 0.0))) {
                throw std::invalid_argument("invalid solvent dielectric parameters");
            }
            if (parameters_.electrolyte_model
                == ELECTROLYTE_MODEL_FIGIEL_SSM_DS) {
                if (!(parameters_.dielectric_ion_suppression_coefficient > 0.0)
                    || !(parameters_.ionic_region_relative_permittivity > 0.0)
                    || !(parameters_.solvation_factor[component] > 0.0)) {
                    throw std::invalid_argument("invalid Figiel electrolyte parameters");
                }
                if (parameters_.charge_number[component] != 0
                    && !(parameters_.born_diameter_angstrom[component] > 0.0)) {
                    throw std::invalid_argument("invalid fitted Born diameter");
                }
            } else if (parameters_.electrolyte_model
                           == ELECTROLYTE_MODEL_ASCANI_ORIGINAL_BORN
                       && parameters_.charge_number[component] != 0
                       && !(parameters_.relative_permittivity[component] > 0.0)) {
                throw std::invalid_argument("invalid ionic dielectric parameter");
            } else if (parameters_.electrolyte_model
                           == ELECTROLYTE_MODEL_HELD_2008_DEBYE_HUCKEL
                       && ((parameters_.charge_number[component] != 0
                            && parameters_.relative_permittivity[component] != 0.0)
                           || parameters_.born_diameter_angstrom[component] != 0.0
                           || parameters_.solvation_factor[component] != 0.0)) {
                throw std::invalid_argument("invalid Held electrolyte parameters");
            }
        }
        if (parameters_.electrolyte_model
                == ELECTROLYTE_MODEL_HELD_2008_DEBYE_HUCKEL
            && (parameters_.dielectric_ion_suppression_coefficient != 0.0
                || parameters_.ionic_region_relative_permittivity != 0.0
                || std::any_of(
                    parameters_.l_ij.begin(),
                    parameters_.l_ij.end(),
                    [](double value) { return value != 0.0; }
                ))) {
            throw std::invalid_argument("invalid Held electrolyte parameters");
        }
    } else if (parameters_.electrolyte_model != ELECTROLYTE_MODEL_NONE) {
        throw std::invalid_argument("neutral model has electrolyte selection");
    }
}

bool NativeModel::supports_native_sdk_v1() const {
    return detail::supports_pure_phase_block(parameters_);
}

bool NativeModel::supports_mixture_native_sdk_v1() const {
    return detail::supports_mixture_phase_block(parameters_);
}

bool NativeModel::supports_electrolyte_native_sdk_v1() const {
    return detail::supports_electrolyte_phase_block(parameters_);
}

std::size_t NativeModel::component_count() const {
    return parameters_.segment_count.size();
}

const std::vector<int>& NativeModel::charge_numbers() const {
    return parameters_.charge_number;
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

NativeMixturePhaseBlock NativeModel::evaluate_mixture_phase(
    double temperature,
    const std::vector<double>& amounts,
    double volume
) const {
    return detail::mixture_phase_block(temperature, amounts, volume, parameters_);
}

NativeMixturePhaseBlock NativeModel::evaluate_mixture_phase_kij(
    double temperature,
    const std::vector<double>& amounts,
    double volume,
    double k_ij
) const {
    return detail::mixture_phase_block_kij(
        temperature, amounts, volume, k_ij, parameters_
    );
}

NativeMixturePhaseBlock NativeModel::evaluate_electrolyte_phase(
    double temperature,
    const std::vector<double>& amounts,
    double volume,
    double ion_mole_fraction_max
) const {
    return detail::electrolyte_phase_block(
        temperature, amounts, volume, ion_mole_fraction_max, parameters_
    );
}

NativePackingFractionBlock NativeModel::evaluate_packing_fraction(
    double temperature,
    const std::vector<double>& amounts,
    double volume,
    double ion_mole_fraction_max
) const {
    return detail::packing_fraction_block(
        temperature, amounts, volume, ion_mole_fraction_max, parameters_
    );
}

std::vector<double> NativeModel::reference_composition(
    const std::vector<double>& solvent_fractions,
    const std::array<std::size_t, 2>& ion_indices,
    const std::array<int, 2>& stoichiometry,
    double molality
) const {
    return detail::formula_composition(
        solvent_fractions,
        parameters_.molar_mass_kg_per_mol,
        ion_indices,
        stoichiometry,
        molality
    );
}

NativeReferenceSequence NativeModel::electrolyte_reference_sequence(
    double temperature,
    double pressure,
    const std::vector<double>& solvent_fractions,
    const std::array<std::size_t, 2>& ion_indices,
    const std::array<int, 2>& stoichiometry,
    const std::string& phase
) const {
    const auto result = detail::electrolyte_reference_sequence(
        solvent_fractions,
        parameters_.molar_mass_kg_per_mol,
        ion_indices,
        stoichiometry,
        [this, temperature, pressure, &ion_indices, &phase](
            const std::vector<double>& composition
        ) {
            const auto state = evaluate_pressure(
                temperature, pressure, composition, phase
            );
            if (!state.log_fugacity_coefficient.has_value()) {
                throw std::runtime_error("reference fugacity is unavailable");
            }
            return detail::ReferenceObservation{
                {
                    (*state.log_fugacity_coefficient)[ion_indices[0]],
                    (*state.log_fugacity_coefficient)[ion_indices[1]],
                },
                std::nullopt,
            };
        }
    );
    return {
        result.log_fugacity,
        result.molality,
        result.convergence_error,
    };
}

NativeIonSolvationBorn NativeModel::evaluate_ion_solvation_born(
    double temperature,
    double pressure,
    double born_diameter_angstrom
) const {
    if (!std::isfinite(temperature) || !std::isfinite(pressure)
        || !std::isfinite(born_diameter_angstrom) || temperature <= 0.0
        || pressure <= 0.0 || born_diameter_angstrom <= 0.0) {
        throw std::invalid_argument("invalid active Born reference inputs");
    }
    if (parameters_.charge_number.size() != 3
        || parameters_.charge_number[0] != 0
        || parameters_.charge_number[1] == 0
        || parameters_.charge_number[2] == 0
        || parameters_.charge_number[1] + parameters_.charge_number[2] != 0) {
        throw std::logic_error("unsupported active Born reference model");
    }
    NativeParameters active_parameters = parameters_;
    active_parameters.born_diameter_angstrom[1] = born_diameter_angstrom;
    const NativeModel active_model(active_parameters);
    const auto result = detail::electrolyte_reference_sequence(
        std::vector<double>{1.0, 0.0, 0.0},
        active_parameters.molar_mass_kg_per_mol,
        std::array<std::size_t, 2>{1, 2},
        std::array<int, 2>{1, 1},
        [&active_model, &active_parameters, temperature, pressure](
            const std::vector<double>& composition
        ) {
            return active_born_observation(
                active_model,
                active_parameters,
                temperature,
                pressure,
                composition,
                1
            );
        }
    );
    return {
        detail::GAS_CONSTANT * temperature * result.log_fugacity[0],
        (*result.derivative)[0],
        result.molality,
        result.convergence_error,
    };
}

std::array<double, 2> NativeModel::molar_volume_bounds(
    double temperature,
    const std::vector<double>& mole_fractions,
    double packing_fraction_min,
    double packing_fraction_max,
    double ion_mole_fraction_max
) const {
    const std::size_t components = parameters_.segment_count.size();
    if (!std::isfinite(temperature) || temperature <= 0.0
        || mole_fractions.size() != components
        || !std::isfinite(packing_fraction_min)
        || !std::isfinite(packing_fraction_max)
        || packing_fraction_min <= 0.0
        || packing_fraction_min >= packing_fraction_max) {
        throw std::invalid_argument("invalid molar-volume bound inputs");
    }
    if (packing_fraction_max >= detail::close_packing_fraction()) {
        throw std::invalid_argument(
            "packing-fraction upper bound is outside the model domain"
        );
    }

    double total_fraction = 0.0;
    double charge_fraction = 0.0;
    double ion_fraction = 0.0;
    for (std::size_t component = 0; component < components; ++component) {
        const double fraction = mole_fractions[component];
        if (!std::isfinite(fraction) || fraction < 0.0) {
            throw std::invalid_argument("invalid model composition");
        }
        total_fraction += fraction;
        const int charge = parameters_.charge_number[component];
        charge_fraction += fraction * static_cast<double>(charge);
        if (charge != 0) {
            ion_fraction += fraction;
        }
    }
    if (std::abs(total_fraction - 1.0) > 1.0e-12
        || std::abs(charge_fraction) > 1.0e-12
        || (std::isfinite(ion_mole_fraction_max)
            && ion_fraction > ion_mole_fraction_max + 1.0e-12)) {
        throw std::invalid_argument("composition is outside the model domain");
    }

    const double coefficient = detail::packing_fraction_coefficient(
        temperature, mole_fractions, parameters_
    );
    if (!(coefficient > 0.0) || !std::isfinite(coefficient)) {
        throw std::invalid_argument("invalid model packing coefficient");
    }
    const std::array<double, 2> bounds = {
        coefficient / packing_fraction_max,
        coefficient / packing_fraction_min,
    };
    if (!(bounds[0] > 0.0) || !(bounds[0] < bounds[1])
        || !std::isfinite(bounds[1])) {
        throw std::invalid_argument("invalid model molar-volume bounds");
    }
    return bounds;
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
        diameter,
        diameter
    );
    const double born_ssm_ds = detail::born_contribution(
        temperature,
        mole_fractions,
        parameters_,
        electrolyte.bulk_relative_permittivity,
        sigma
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
    const double packing_coefficient = detail::packing_fraction_coefficient(
        temperature, mole_fractions, parameters_
    );
    if (!(packing_coefficient > 0.0) || !std::isfinite(packing_coefficient)) {
        throw std::invalid_argument("invalid pressure-state packing coefficient");
    }
    const auto root = detail::solve_density(
        pressure,
        pressure / (detail::GAS_CONSTANT * temperature),
        detail::close_packing_fraction() / packing_coefficient,
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
