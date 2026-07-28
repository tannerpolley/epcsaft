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

// EqExclude: input and domain validation
bool positive_finite(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value) && value > 0.0;
    });
}

// EqExclude: input and domain validation
bool finite(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
    });
}

double log_fugacity_density(
    const NativeResult& state,
    const detail::ThermodynamicDerivatives& derivatives,
    const std::vector<double>& composition,
    std::size_t component,
    std::size_t coordinates
) {
    // EqID: log_fugacity_density_derivative
    auto hessian = [&derivatives, coordinates](
        std::size_t row, std::size_t column
    ) {
        return derivatives.hessian[row * coordinates + column];
    };
    double weighted_composition_density = 0.0;
    for (std::size_t mixture_component = 0;
         mixture_component < composition.size();
         ++mixture_component) {
        weighted_composition_density += composition[mixture_component]
            * hessian(2 + mixture_component, 1);
    }
    const double compressibility_density = derivatives.density
        + state.molar_density * derivatives.density_density;
    const double chemical_potential_density = 2.0 * derivatives.density
        + state.molar_density * derivatives.density_density
        + hessian(2 + component, 1)
        - weighted_composition_density;
    return chemical_potential_density
        - compressibility_density / state.compressibility_factor;
}

double fixed_pressure_log_fugacity(
    const NativeResult& state,
    double pressure,
    std::size_t component,
    double fugacity_density
) {
    // EqID: fixed_pressure_fugacity_correction
    const double density_correction = -(state.pressure - pressure)
        / state.pressure_density_derivative;
    return (*state.log_fugacity_coefficient)[component]
        + fugacity_density * density_correction;
}

enum class ActiveSolvationParameter {
    BornDiameter,
    IonicRegionRelativePermittivity,
    SolventRelativePermittivity,
};

// EqID: fixed_pressure_total_derivative
detail::ReferenceObservation active_solvation_parameter_observation(
    const NativeModel& model,
    const NativeParameters& parameters,
    double temperature,
    double pressure,
    const std::vector<double>& composition,
    std::size_t active_component,
    ActiveSolvationParameter active_parameter
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
    const auto derivatives = [&] {
        if (active_parameter == ActiveSolvationParameter::BornDiameter) {
            return detail::thermodynamic_derivatives_active_born(
                temperature,
                state.molar_density,
                composition,
                parameters,
                site_fractions,
                active_component,
                parameters.born_diameter_angstrom[active_component]
            );
        }
        if (active_parameter
            == ActiveSolvationParameter::IonicRegionRelativePermittivity) {
            return detail::
                thermodynamic_derivatives_active_ionic_region_permittivity(
                    temperature,
                    state.molar_density,
                    composition,
                    parameters,
                    site_fractions,
                    parameters.ionic_region_relative_permittivity
                );
        }
        return detail::thermodynamic_derivatives_active_relative_permittivity(
            temperature,
            state.molar_density,
            composition,
            parameters,
            site_fractions,
            active_component,
            parameters.relative_permittivity[active_component]
        );
    }();
    const std::optional<double>& active_derivative =
        active_parameter == ActiveSolvationParameter::BornDiameter
        ? derivatives.active_born
        : active_parameter
                == ActiveSolvationParameter::IonicRegionRelativePermittivity
            ? derivatives.active_ionic_permittivity
            : derivatives.active_relative_permittivity;
    if (!active_derivative.has_value()) {
        throw std::runtime_error("active solvation derivative is unavailable");
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
    std::vector<double> total_derivative(2);
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
        const double chemical_potential_diameter = *active_derivative
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

// EqID: infinite_dilution_solvation_energy
// EqID: fixed_pressure_total_derivative
NativeIonSolvationIonicPermittivity active_solvation_parameter_reference(
    NativeParameters active_parameters,
    double temperature,
    double pressure,
    std::size_t active_component,
    ActiveSolvationParameter active_parameter
) {
    const NativeModel active_model(active_parameters);
    const auto result = detail::electrolyte_reference_sequence(
        std::vector<double>{1.0, 0.0, 0.0},
        active_parameters.molar_mass_kg_per_mol,
        std::array<std::size_t, 2>{1, 2},
        std::array<int, 2>{1, 1},
        [&active_model,
         &active_parameters,
         temperature,
         pressure,
         active_component,
         active_parameter](const std::vector<double>& composition) {
            return active_solvation_parameter_observation(
                active_model,
                active_parameters,
                temperature,
                pressure,
                composition,
                active_component,
                active_parameter
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

// EqID: fixed_pressure_total_derivative
detail::ReferenceObservation active_aqueous_kij_observation(
    const NativeModel& model,
    const NativeParameters& parameters,
    double temperature,
    double pressure,
    const std::vector<double>& composition,
    const std::array<double, 3>& k_ij
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
    const auto derivatives = detail::thermodynamic_derivatives_active_aqueous_kij(
        temperature,
        state.molar_density,
        composition,
        parameters,
        site_fractions,
        k_ij
    );
    if (derivatives.active_aqueous_k_ij.size() != 3) {
        throw std::runtime_error("active aqueous k_ij derivative is unavailable");
    }
    const std::size_t components = composition.size();
    const std::size_t coordinates = 2 + components + 3;
    const std::size_t active_offset = 2 + components;
    const auto hessian = [&derivatives, coordinates](
        std::size_t row, std::size_t column
    ) {
        return derivatives.hessian[row * coordinates + column];
    };
    const double compressibility = state.compressibility_factor;
    std::array<double, 2> fugacity_density{};
    std::array<double, 2> corrected_log_fugacity{};
    for (std::size_t ion = 0; ion < 2; ++ion) {
        fugacity_density[ion] = log_fugacity_density(
            state, derivatives, composition, 1 + ion, coordinates
        );
        corrected_log_fugacity[ion] = fixed_pressure_log_fugacity(
            state, pressure, 1 + ion, fugacity_density[ion]
        );
    }
    std::vector<double> total_derivative(6);
    for (std::size_t active = 0; active < 3; ++active) {
        const std::size_t active_coordinate = active_offset + active;
        const double pressure_k = state.molar_density * state.molar_density
            * detail::GAS_CONSTANT * temperature * hessian(1, active_coordinate);
        const double density_k = -pressure_k / state.pressure_density_derivative;
        for (std::size_t ion = 0; ion < 2; ++ion) {
            const std::size_t component = 1 + ion;
            double weighted_composition_k = 0.0;
            for (std::size_t mixture_component = 0;
                 mixture_component < components;
                 ++mixture_component) {
                weighted_composition_k += composition[mixture_component]
                    * hessian(2 + mixture_component, active_coordinate);
            }
            const double chemical_potential_k =
                derivatives.active_aqueous_k_ij[active]
                + state.molar_density * hessian(1, active_coordinate)
                + hessian(2 + component, active_coordinate)
                - weighted_composition_k;
            const double log_fugacity_k = chemical_potential_k
                - state.molar_density * hessian(1, active_coordinate)
                    / compressibility;
            total_derivative[ion * 3 + active] = log_fugacity_k
                + fugacity_density[ion] * density_k;
        }
    }
    return {
        corrected_log_fugacity,
        std::move(total_derivative),
    };
}

// EqID: fixed_pressure_total_derivative
detail::ReferenceObservation active_solvation_factor_observation(
    const NativeModel& model,
    const NativeParameters& parameters,
    double temperature,
    double pressure,
    const std::vector<double>& composition,
    double solvation_factor
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
    const auto derivatives =
        detail::thermodynamic_derivatives_active_solvation_factor(
            temperature,
            state.molar_density,
            composition,
            parameters,
            site_fractions,
            solvation_factor
        );
    if (!derivatives.active_solvation_factor.has_value()) {
        throw std::runtime_error("active solvation-factor derivative is unavailable");
    }
    const std::size_t components = composition.size();
    const std::size_t coordinates = 2 + components + 1;
    const std::size_t active_coordinate = 2 + components;
    const auto hessian = [&derivatives, coordinates](
        std::size_t row, std::size_t column
    ) {
        return derivatives.hessian[row * coordinates + column];
    };
    const double pressure_factor = state.molar_density * state.molar_density
        * detail::GAS_CONSTANT * temperature * hessian(1, active_coordinate);
    const double density_factor = -pressure_factor
        / state.pressure_density_derivative;
    std::array<double, 2> corrected_log_fugacity{};
    std::vector<double> total_derivative(2);
    for (std::size_t ion = 0; ion < 2; ++ion) {
        const std::size_t component = 1 + ion;
        const double fugacity_density = log_fugacity_density(
            state, derivatives, composition, component, coordinates
        );
        corrected_log_fugacity[ion] = fixed_pressure_log_fugacity(
            state, pressure, component, fugacity_density
        );
        double weighted_composition_factor = 0.0;
        for (std::size_t mixture_component = 0;
             mixture_component < components;
             ++mixture_component) {
            weighted_composition_factor += composition[mixture_component]
                * hessian(2 + mixture_component, active_coordinate);
        }
        const double chemical_potential_factor =
            *derivatives.active_solvation_factor
            + state.molar_density * hessian(1, active_coordinate)
            + hessian(2 + component, active_coordinate)
            - weighted_composition_factor;
        const double log_fugacity_factor = chemical_potential_factor
            - state.molar_density * hessian(1, active_coordinate)
                / state.compressibility_factor;
        total_derivative[ion] = log_fugacity_factor
            + fugacity_density * density_factor;
    }
    return {corrected_log_fugacity, std::move(total_derivative)};
}

// EqExclude: transport and buffer mechanics
std::optional<double> single_solvent_log_fugacity(
    const NativeModel& model,
    double temperature,
    double pressure,
    const std::vector<double>& solvent_fractions,
    const std::string& phase
) {
    std::optional<std::size_t> solvent_index;
    for (std::size_t component = 0; component < solvent_fractions.size(); ++component) {
        if (solvent_fractions[component] <= 0.0) {
            continue;
        }
        if (solvent_index.has_value()) {
            return std::nullopt;
        }
        solvent_index = component;
    }
    if (!solvent_index.has_value()) {
        return std::nullopt;
    }
    const auto state = model.evaluate_pressure(
        temperature, pressure, solvent_fractions, phase
    );
    if (!state.log_fugacity_coefficient.has_value()) {
        throw std::runtime_error("pure-solvent reference fugacity is unavailable");
    }
    return (*state.log_fugacity_coefficient)[*solvent_index];
}

}  // namespace

// EqExclude: input and domain validation
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
        parameters_.solvation_factor.size() != components ||
        parameters_.dielectric_volume_m3_per_mol.size() != components ||
        parameters_.dielectric_temperature_k.size() != components ||
        parameters_.ion_specific_suppression_coefficient.size() != components ||
        parameters_.dipole_scaling_debye2.size() != components ||
        parameters_.polarizability_scaling_angstrom3.size() != components ||
        parameters_.correlation_integral_parameter.size() != components) {
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
        !std::isfinite(parameters_.ion_fraction_suppression_coefficient) ||
        !finite(parameters_.born_diameter_angstrom) ||
        !finite(parameters_.solvation_factor) ||
        !finite(parameters_.dielectric_volume_m3_per_mol) ||
        !finite(parameters_.dielectric_temperature_k) ||
        !finite(parameters_.ion_specific_suppression_coefficient) ||
        !finite(parameters_.dipole_scaling_debye2) ||
        !finite(parameters_.polarizability_scaling_angstrom3) ||
        !finite(parameters_.correlation_integral_parameter) ||
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
        if (parameters_.electrolyte_model != ELECTROLYTE_MODEL_SSM_DS
            && parameters_.electrolyte_model
                != ELECTROLYTE_MODEL_ORIGINAL_BORN
            && parameters_.electrolyte_model
                != ELECTROLYTE_MODEL_FULLY_DISSOCIATED_DEBYE_HUCKEL) {
            throw std::invalid_argument("invalid electrolyte model selection");
        }
        const bool ssm_ds = parameters_.electrolyte_model
            == ELECTROLYTE_MODEL_SSM_DS;
        const bool original_born = parameters_.electrolyte_model
            == ELECTROLYTE_MODEL_ORIGINAL_BORN;
        const bool fully_dissociated = parameters_.electrolyte_model
            == ELECTROLYTE_MODEL_FULLY_DISSOCIATED_DEBYE_HUCKEL;
        const bool valid_ssm_permittivity = parameters_.relative_permittivity_formulation
                == RELATIVE_PERMITTIVITY_FORMULATION_ION_FRACTION_SUPPRESSION
            || parameters_.relative_permittivity_formulation
                == RELATIVE_PERMITTIVITY_FORMULATION_DENSITY_VOLUME_TEMPERATURE
            || parameters_.relative_permittivity_formulation
                == RELATIVE_PERMITTIVITY_FORMULATION_ION_SPECIFIC_SUPPRESSION
            || parameters_.relative_permittivity_formulation
                == RELATIVE_PERMITTIVITY_FORMULATION_DIPOLAR_DENSITY_PERTURBATION;
        if ((ssm_ds && !valid_ssm_permittivity)
            || (original_born && parameters_.relative_permittivity_formulation
                != RELATIVE_PERMITTIVITY_FORMULATION_COMPONENT_MIXING)
            || (fully_dissociated && parameters_.relative_permittivity_formulation
                != RELATIVE_PERMITTIVITY_FORMULATION_SOLVENT_ONLY)) {
            throw std::invalid_argument("incompatible electrolyte and permittivity models");
        }
        std::size_t solvent_count = 0;
        for (std::size_t component = 0; component < components; ++component) {
            const bool solvent = parameters_.charge_number[component] == 0;
            solvent_count += static_cast<std::size_t>(solvent);
            const bool needs_tabulated_permittivity = parameters_.relative_permittivity_formulation
                    == RELATIVE_PERMITTIVITY_FORMULATION_ION_FRACTION_SUPPRESSION
                || parameters_.relative_permittivity_formulation
                    == RELATIVE_PERMITTIVITY_FORMULATION_COMPONENT_MIXING
                || parameters_.relative_permittivity_formulation
                    == RELATIVE_PERMITTIVITY_FORMULATION_ION_SPECIFIC_SUPPRESSION
                || parameters_.relative_permittivity_formulation
                    == RELATIVE_PERMITTIVITY_FORMULATION_SOLVENT_ONLY;
            if (solvent
                && (!(parameters_.molar_mass_kg_per_mol[component] > 0.0)
                    || (needs_tabulated_permittivity
                        && !(parameters_.relative_permittivity[component] > 0.0)))) {
                throw std::invalid_argument("invalid solvent dielectric parameters");
            }
            if (ssm_ds) {
                if (!(parameters_.ionic_region_relative_permittivity > 0.0)
                    || !(parameters_.solvation_factor[component] > 0.0)) {
                    throw std::invalid_argument("invalid SSM+DS electrolyte parameters");
                }
                if (!solvent
                    && !(parameters_.born_diameter_angstrom[component] > 0.0)) {
                    throw std::invalid_argument("invalid fitted Born diameter");
                }
            } else if (original_born && !solvent
                       && !(parameters_.relative_permittivity[component] > 0.0)) {
                throw std::invalid_argument("invalid ionic dielectric parameter");
            } else if (parameters_.electrolyte_model
                           == ELECTROLYTE_MODEL_FULLY_DISSOCIATED_DEBYE_HUCKEL
                       && ((parameters_.charge_number[component] != 0
                            && parameters_.relative_permittivity[component] != 0.0)
                           || parameters_.born_diameter_angstrom[component] != 0.0
                           || parameters_.solvation_factor[component] != 0.0)) {
                throw std::invalid_argument("invalid fully dissociated electrolyte parameters");
            }
            if (parameters_.relative_permittivity_formulation
                    == RELATIVE_PERMITTIVITY_FORMULATION_DENSITY_VOLUME_TEMPERATURE
                && solvent
                && (!(parameters_.dielectric_volume_m3_per_mol[component] > 0.0)
                    || !(parameters_.dielectric_temperature_k[component] > 0.0))) {
                throw std::invalid_argument("invalid density formulation solvent parameters");
            }
            if (parameters_.relative_permittivity_formulation
                    == RELATIVE_PERMITTIVITY_FORMULATION_ION_SPECIFIC_SUPPRESSION
                && !solvent
                && !(parameters_.ion_specific_suppression_coefficient[component] > 0.0)) {
                throw std::invalid_argument("invalid ion-specific suppression parameter");
            }
            if (parameters_.relative_permittivity_formulation
                == RELATIVE_PERMITTIVITY_FORMULATION_DIPOLAR_DENSITY_PERTURBATION) {
                if (parameters_.dipole_scaling_debye2[component] < 0.0
                    || parameters_.polarizability_scaling_angstrom3[component] < 0.0
                    || !(parameters_.correlation_integral_parameter[component] > 0.0)) {
                    throw std::invalid_argument("invalid dipolar perturbation parameters");
                }
            }
        }
        if (parameters_.relative_permittivity_formulation == RELATIVE_PERMITTIVITY_FORMULATION_ION_FRACTION_SUPPRESSION
            && !(parameters_.ion_fraction_suppression_coefficient > 0.0)) {
            throw std::invalid_argument("invalid ion-fraction suppression parameter");
        }
        if (parameters_.relative_permittivity_formulation == RELATIVE_PERMITTIVITY_FORMULATION_ION_SPECIFIC_SUPPRESSION
            && solvent_count != 1) {
            throw std::invalid_argument("mixed-solvent volume data are required");
        }
        if (fully_dissociated
            && (parameters_.ion_fraction_suppression_coefficient != 0.0
                || parameters_.ionic_region_relative_permittivity != 0.0
                || std::any_of(
                    parameters_.l_ij.begin(),
                    parameters_.l_ij.end(),
                    [](double value) { return value != 0.0; }
                ))) {
            throw std::invalid_argument("invalid fully dissociated electrolyte parameters");
        }
    } else if (parameters_.electrolyte_model != ELECTROLYTE_MODEL_NONE
               || parameters_.relative_permittivity_formulation != RELATIVE_PERMITTIVITY_FORMULATION_NONE) {
        throw std::invalid_argument("neutral model has electrolyte selection");
    }
}

// EqExclude: input and domain validation
bool NativeModel::supports_native_sdk_v1() const {
    return detail::supports_pure_phase_block(parameters_);
}

// EqExclude: input and domain validation
bool NativeModel::supports_parameterized_pure_native_sdk_v1() const {
    return detail::supports_pure_phase_block(parameters_)
        && parameters_.sigma_term_component.empty()
        && parameters_.sigma_term_amplitude.empty()
        && parameters_.sigma_term_exponent.empty();
}

// EqExclude: input and domain validation
bool NativeModel::supports_associating_pure_parameter_native_sdk_v1() const {
    return detail::supports_associating_pure_phase_block(parameters_);
}

// EqExclude: input and domain validation
bool NativeModel::supports_mixture_native_sdk_v1() const {
    return detail::supports_mixture_phase_block(parameters_);
}

// EqExclude: input and domain validation
bool NativeModel::supports_electrolyte_native_sdk_v1() const {
    return detail::supports_electrolyte_phase_block(parameters_);
}

// EqExclude: serialization
std::size_t NativeModel::component_count() const {
    return parameters_.segment_count.size();
}

// EqExclude: serialization
const std::vector<int>& NativeModel::charge_numbers() const {
    return parameters_.charge_number;
}

// EqExclude: serialization
const std::vector<double>& NativeModel::molar_masses_kg_per_mol() const {
    return parameters_.molar_mass_kg_per_mol;
}

// EqExclude: transport and buffer mechanics
NativePurePhaseBlock NativeModel::evaluate_pure_phase(
    double temperature,
    double amount,
    double volume
) const {
    return detail::pure_phase_block(temperature, amount, volume, parameters_);
}

// EqExclude: transport and buffer mechanics
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

// EqExclude: transport and buffer mechanics
NativeAssociatingParameterizedPurePhaseBlock
NativeModel::evaluate_diameter_basis_associating_pure_phase_parameters(
    double temperature,
    double amount,
    double volume,
    double segment_count,
    double segment_diameter_angstrom,
    double dispersion_energy_over_k,
    double association_energy_over_k,
    double association_volume
) const {
    return detail::associating_parameterized_pure_phase_block(
        temperature,
        amount,
        volume,
        segment_count,
        segment_diameter_angstrom,
        dispersion_energy_over_k,
        association_energy_over_k,
        association_volume,
        parameters_
    );
}

// EqExclude: transport and buffer mechanics
NativeMixturePhaseBlock NativeModel::evaluate_pure_phase_parameter(
    double temperature,
    double amount,
    double volume,
    ActivePureParameter family,
    double parameter_value
) const {
    if (family == ActivePureParameter::AssociationEnergyOverK
        || family == ActivePureParameter::AssociationVolume) {
        return detail::associating_pure_phase_parameter_block(
            temperature,
            amount,
            volume,
            parameter_value,
            family == ActivePureParameter::AssociationEnergyOverK,
            parameters_
        );
    }
    if (!supports_parameterized_pure_native_sdk_v1()) {
        throw std::logic_error(
            "active pure parameters require one neutral nonassociating "
            "component with constant segment diameter"
        );
    }
    double segment_count = parameters_.segment_count.at(0);
    double segment_diameter = parameters_.sigma_constant.at(0);
    double dispersion_energy = parameters_.dispersion_energy_over_k.at(0);
    std::size_t active_coordinate = 0;
    switch (family) {
        case ActivePureParameter::SegmentCount:
            segment_count = parameter_value;
            active_coordinate = 2;
            break;
        case ActivePureParameter::SegmentDiameter:
            segment_diameter = parameter_value;
            active_coordinate = 3;
            break;
        case ActivePureParameter::DispersionEnergyOverK:
            dispersion_energy = parameter_value;
            active_coordinate = 4;
            break;
        case ActivePureParameter::AssociationEnergyOverK:
        case ActivePureParameter::AssociationVolume:
            throw std::logic_error("association parameter dispatch failed");
        default:
            throw std::invalid_argument("unknown active pure parameter");
    }
    const auto full = evaluate_pure_phase_parameters(
        temperature,
        amount,
        volume,
        segment_count,
        segment_diameter,
        dispersion_energy
    );
    const std::array<std::size_t, 3> selected{0, 1, active_coordinate};
    std::vector<double> gradient(3);
    std::vector<double> hessian(9);
    for (std::size_t row = 0; row < selected.size(); ++row) {
        gradient[row] = full.gradient[selected[row]];
        for (std::size_t column = 0; column < selected.size(); ++column) {
            hessian[row * selected.size() + column] =
                full.hessian[selected[row] * full.gradient.size()
                             + selected[column]];
        }
    }
    return {
        full.helmholtz_over_rt_reference_amount,
        std::move(gradient),
        std::move(hessian),
        full.pressure_pa,
    };
}

// EqExclude: transport and buffer mechanics
NativeMixturePhaseBlock NativeModel::evaluate_mixture_phase(
    double temperature,
    const std::vector<double>& amounts,
    double volume
) const {
    return detail::mixture_phase_block(temperature, amounts, volume, parameters_);
}

// EqExclude: transport and buffer mechanics
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

// EqExclude: transport and buffer mechanics
NativeMixturePhaseBlock NativeModel::evaluate_mixture_phase_lij(
    double temperature,
    const std::vector<double>& amounts,
    double volume,
    double l_ij
) const {
    return detail::mixture_phase_block_lij(
        temperature, amounts, volume, l_ij, parameters_
    );
}

// EqExclude: transport and buffer mechanics
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

// EqExclude: transport and buffer mechanics
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

// EqID: formula_unit_composition
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

// EqID: reference_molality_sequence
NativeReferenceSequence NativeModel::electrolyte_reference_sequence(
    double temperature,
    double pressure,
    const std::vector<double>& solvent_fractions,
    const std::array<std::size_t, 2>& ion_indices,
    const std::array<int, 2>& stoichiometry,
    const std::string& phase,
    NativeEvaluationBudget* budget
) const {
    auto checkpoint = [budget] {
        if (budget != nullptr) {
            budget->checkpoint();
        }
    };
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
        },
        checkpoint
    );
    const auto pure_solvent_log_fugacity = single_solvent_log_fugacity(
        *this, temperature, pressure, solvent_fractions, phase
    );
    return {
        result.log_fugacity,
        pure_solvent_log_fugacity,
        result.molality,
        result.convergence_error,
    };
}

// EqID: infinite_dilution_solvation_energy
// EqID: fixed_pressure_total_derivative
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
    const auto result = active_solvation_parameter_reference(
        std::move(active_parameters),
        temperature,
        pressure,
        1,
        ActiveSolvationParameter::BornDiameter
    );
    return {
        result.solvation_gibbs_j_per_mol,
        result.derivative_j_per_mol,
        result.reference_molality_mol_per_kg,
        result.reference_convergence_error,
    };
}

// EqID: infinite_dilution_solvation_energy
// EqID: fixed_pressure_total_derivative
NativeIonSolvationIonicPermittivity
NativeModel::evaluate_ion_solvation_ionic_permittivity(
    double temperature,
    double pressure,
    double ionic_region_relative_permittivity
) const {
    if (!std::isfinite(temperature) || !std::isfinite(pressure)
        || !std::isfinite(ionic_region_relative_permittivity)
        || temperature <= 0.0 || pressure <= 0.0
        || ionic_region_relative_permittivity <= 1.0) {
        throw std::invalid_argument(
            "invalid active ionic-region permittivity inputs"
        );
    }
    if (parameters_.charge_number.size() != 3
        || parameters_.charge_number[0] != 0
        || parameters_.charge_number[1] == 0
        || parameters_.charge_number[2] == 0
        || parameters_.charge_number[1] + parameters_.charge_number[2] != 0) {
        throw std::logic_error(
            "unsupported active ionic-region permittivity model"
        );
    }
    NativeParameters active_parameters = parameters_;
    active_parameters.ionic_region_relative_permittivity =
        ionic_region_relative_permittivity;
    return active_solvation_parameter_reference(
        std::move(active_parameters),
        temperature,
        pressure,
        1,
        ActiveSolvationParameter::IonicRegionRelativePermittivity
    );
}

// EqID: infinite_dilution_solvation_energy
// EqID: fixed_pressure_total_derivative
NativeIonSolvationIonicPermittivity
NativeModel::evaluate_ion_solvation_solvent_permittivity(
    double temperature,
    double pressure,
    double solvent_relative_permittivity
) const {
    if (!std::isfinite(temperature) || !std::isfinite(pressure)
        || !std::isfinite(solvent_relative_permittivity)
        || temperature <= 0.0 || pressure <= 0.0
        || solvent_relative_permittivity <= 1.0) {
        throw std::invalid_argument(
            "invalid active solvent relative-permittivity inputs"
        );
    }
    if (parameters_.charge_number.size() != 3
        || parameters_.charge_number[0] != 0
        || parameters_.charge_number[1] == 0
        || parameters_.charge_number[2] == 0
        || parameters_.charge_number[1] + parameters_.charge_number[2] != 0
        || parameters_.relative_permittivity.size() != 3) {
        throw std::logic_error(
            "unsupported active solvent relative-permittivity model"
        );
    }
    NativeParameters active_parameters = parameters_;
    active_parameters.relative_permittivity[0] =
        solvent_relative_permittivity;
    return active_solvation_parameter_reference(
        std::move(active_parameters),
        temperature,
        pressure,
        0,
        ActiveSolvationParameter::SolventRelativePermittivity
    );
}

// EqID: infinite_dilution_solvation_energy
// EqID: fixed_pressure_total_derivative
NativeIonSolvationKij NativeModel::evaluate_ion_solvation_kij(
    double temperature,
    double pressure,
    std::size_t active_ion_component,
    const std::array<double, 3>& k_ij
) const {
    if (!std::isfinite(temperature) || !std::isfinite(pressure)
        || temperature <= 0.0 || pressure <= 0.0
        || (active_ion_component != 1 && active_ion_component != 2)
        || !std::all_of(k_ij.begin(), k_ij.end(), [](double value) {
            return std::isfinite(value);
        })) {
        throw std::invalid_argument("invalid active ion-solvation k_ij inputs");
    }
    if (parameters_.charge_number.size() != 3
        || parameters_.charge_number[0] != 0
        || parameters_.charge_number[1] == 0
        || parameters_.charge_number[2] == 0
        || parameters_.charge_number[1] + parameters_.charge_number[2] != 0) {
        throw std::logic_error("unsupported active ion-solvation k_ij model");
    }
    NativeParameters active_parameters = parameters_;
    active_parameters.k_ij[1] = active_parameters.k_ij[3] = k_ij[0];
    active_parameters.k_ij[2] = active_parameters.k_ij[6] = k_ij[1];
    active_parameters.k_ij[5] = active_parameters.k_ij[7] = k_ij[2];
    const NativeModel active_model(active_parameters);
    const auto result = detail::electrolyte_reference_sequence(
        std::vector<double>{1.0, 0.0, 0.0},
        active_parameters.molar_mass_kg_per_mol,
        std::array<std::size_t, 2>{1, 2},
        std::array<int, 2>{1, 1},
        [&active_model, &active_parameters, temperature, pressure, &k_ij](
            const std::vector<double>& composition
        ) {
            return active_aqueous_kij_observation(
                active_model,
                active_parameters,
                temperature,
                pressure,
                composition,
                k_ij
            );
        }
    );
    if (!result.derivative.has_value() || result.derivative->size() != 6) {
        throw std::runtime_error(
            "active ion-solvation k_ij derivative is unavailable"
        );
    }
    const std::size_t ion = active_ion_component - 1;
    std::array<double, 3> derivative{};
    for (std::size_t active = 0; active < 3; ++active) {
        derivative[active] = detail::GAS_CONSTANT * temperature
            * (*result.derivative)[ion * 3 + active];
    }
    return {
        detail::GAS_CONSTANT * temperature * result.log_fugacity[ion],
        derivative,
        result.molality,
        result.convergence_error,
    };
}

// EqExclude: transport and buffer mechanics
NativeAqueousMiacKij NativeModel::evaluate_aqueous_miac_kij(
    double temperature,
    double pressure,
    double molality,
    const std::array<double, 3>& k_ij
) const {
    return evaluate_aqueous_miac_kij_batch(
        temperature, pressure, std::vector<double>{molality}, k_ij
    ).front();
}

std::vector<NativeAqueousMiacKij>
// EqID: mean_ionic_activity_coefficient
// EqID: molality_activity_conversion
// EqID: fixed_pressure_total_derivative
NativeModel::evaluate_aqueous_miac_kij_batch(
    double temperature,
    double pressure,
    const std::vector<double>& molalities,
    const std::array<double, 3>& k_ij,
    NativeEvaluationBudget* budget
) const {
    if (!std::isfinite(temperature) || !std::isfinite(pressure)
        || temperature <= 0.0 || pressure <= 0.0
        || molalities.empty() || !positive_finite(molalities)
        || !std::all_of(k_ij.begin(), k_ij.end(), [](double value) {
            return std::isfinite(value);
        })) {
        throw std::invalid_argument("invalid active aqueous k_ij inputs");
    }
    if (parameters_.charge_number.size() != 3
        || parameters_.charge_number != std::vector<int>({0, 1, -1})) {
        throw std::logic_error("unsupported active aqueous k_ij model");
    }
    NativeParameters active_parameters = parameters_;
    active_parameters.k_ij[1] = active_parameters.k_ij[3] = k_ij[0];
    active_parameters.k_ij[2] = active_parameters.k_ij[6] = k_ij[1];
    active_parameters.k_ij[5] = active_parameters.k_ij[7] = k_ij[2];
    const NativeModel active_model(active_parameters);
    auto checkpoint = [budget] {
        if (budget != nullptr) {
            budget->checkpoint();
        }
    };
    const auto reference = detail::electrolyte_reference_sequence(
        std::vector<double>{1.0, 0.0, 0.0},
        active_parameters.molar_mass_kg_per_mol,
        std::array<std::size_t, 2>{1, 2},
        std::array<int, 2>{1, 1},
        [&active_model, &active_parameters, temperature, pressure, &k_ij](
            const std::vector<double>& reference_composition
        ) {
            return active_aqueous_kij_observation(
                active_model,
                active_parameters,
                temperature,
                pressure,
                reference_composition,
                k_ij
            );
        },
        checkpoint
    );
    if (!reference.derivative.has_value() || reference.derivative->size() != 6
        || !reference.derivative_convergence_error.has_value()) {
        throw std::runtime_error("active aqueous k_ij derivative is unavailable");
    }
    std::vector<NativeAqueousMiacKij> results(molalities.size());
    detail::cppad_parallel_for(molalities.size(), [&](std::size_t index) {
        checkpoint();
        const double molality = molalities[index];
        const auto composition = detail::formula_composition(
            std::vector<double>{1.0, 0.0, 0.0},
            active_parameters.molar_mass_kg_per_mol,
            std::array<std::size_t, 2>{1, 2},
            std::array<int, 2>{1, 1},
            molality
        );
        const auto finite = active_aqueous_kij_observation(
            active_model,
            active_parameters,
            temperature,
            pressure,
            composition,
            k_ij
        );
        if (!finite.derivative.has_value() || finite.derivative->size() != 6) {
            throw std::runtime_error(
                "active aqueous k_ij derivative is unavailable"
            );
        }
        std::array<double, 3> total_derivative{};
        for (std::size_t active = 0; active < 3; ++active) {
            total_derivative[active] = 0.5 * (
                (*finite.derivative)[active] + (*finite.derivative)[3 + active]
                - (*reference.derivative)[active]
                - (*reference.derivative)[3 + active]
            );
        }
        const double conversion = std::log(
            1.0 + 2.0 * active_parameters.molar_mass_kg_per_mol[0] * molality
        );
        results[index] = {
            0.5 * (
                finite.log_fugacity[0] + finite.log_fugacity[1]
                - reference.log_fugacity[0] - reference.log_fugacity[1]
            ) - conversion,
            total_derivative,
            reference.molality,
            reference.convergence_error,
            *reference.derivative_convergence_error,
        };
        const auto& result = results[index];
        if (!std::isfinite(
                result.log_mean_ionic_activity_coefficient_molality
            )
            || !std::all_of(
                result.derivative.begin(),
                result.derivative.end(),
                [](double value) { return std::isfinite(value); }
            )
            || !std::isfinite(result.reference_molality_mol_per_kg)
            || !std::isfinite(result.reference_convergence_error)
            || !std::isfinite(
                result.reference_derivative_convergence_error
            )) {
            throw std::runtime_error(
                "active aqueous k_ij result contains nonfinite values"
            );
        }
    });
    return results;
}

NativeAqueousMiacSolvationFactor
// EqExclude: transport and buffer mechanics
NativeModel::evaluate_aqueous_miac_solvation_factor(
    double temperature,
    double pressure,
    double molality,
    double solvation_factor
) const {
    return evaluate_aqueous_miac_solvation_factor_batch(
        temperature,
        pressure,
        std::vector<double>{molality},
        solvation_factor
    ).front();
}

// EqExclude: transport and buffer mechanics
NativeIonFractionSuppression NativeModel::evaluate_ion_fraction_suppression(
    double total_ion_mole_fraction,
    double suppression_coefficient
) const {
    if (!std::isfinite(total_ion_mole_fraction)
        || !std::isfinite(suppression_coefficient)
        || total_ion_mole_fraction <= 0.0
        || total_ion_mole_fraction >= 1.0
        || suppression_coefficient <= 0.0) {
        throw std::invalid_argument(
            "invalid active ion-fraction suppression inputs"
        );
    }
    if (parameters_.relative_permittivity_formulation
            != RELATIVE_PERMITTIVITY_FORMULATION_ION_FRACTION_SUPPRESSION
        || parameters_.charge_number.size() != 3
        || parameters_.charge_number[0] != 0
        || parameters_.charge_number[1] <= 0
        || parameters_.charge_number[2] >= 0
        || !(parameters_.relative_permittivity[0] > 0.0)) {
        throw std::logic_error(
            "unsupported active ion-fraction suppression model"
        );
    }
    const double salt_free = parameters_.relative_permittivity[0];
    return {
        detail::ion_fraction_suppressed_relative_permittivity(
            salt_free,
            total_ion_mole_fraction,
            suppression_coefficient
        ),
        detail::ion_fraction_suppressed_relative_permittivity_derivative(
            salt_free,
            total_ion_mole_fraction,
            suppression_coefficient
        ),
        salt_free,
    };
}

std::vector<NativeAqueousMiacSolvationFactor>
// EqID: mean_ionic_activity_coefficient
// EqID: molality_activity_conversion
// EqID: fixed_pressure_total_derivative
NativeModel::evaluate_aqueous_miac_solvation_factor_batch(
    double temperature,
    double pressure,
    const std::vector<double>& molalities,
    double solvation_factor,
    NativeEvaluationBudget* budget
) const {
    if (!std::isfinite(temperature) || !std::isfinite(pressure)
        || !std::isfinite(solvation_factor)
        || temperature <= 0.0 || pressure <= 0.0
        || molalities.empty() || !positive_finite(molalities)
        || solvation_factor <= 0.0) {
        throw std::invalid_argument("invalid active solvation-factor inputs");
    }
    if (parameters_.charge_number != std::vector<int>({0, 1, -1})) {
        throw std::logic_error("unsupported active solvation-factor model");
    }
    NativeParameters active_parameters = parameters_;
    active_parameters.solvation_factor[0] = solvation_factor;
    const NativeModel active_model(active_parameters);
    auto checkpoint = [budget] {
        if (budget != nullptr) {
            budget->checkpoint();
        }
    };
    const auto reference = detail::electrolyte_reference_sequence(
        std::vector<double>{1.0, 0.0, 0.0},
        active_parameters.molar_mass_kg_per_mol,
        std::array<std::size_t, 2>{1, 2},
        std::array<int, 2>{1, 1},
        [&active_model, &active_parameters, temperature, pressure, solvation_factor](
            const std::vector<double>& reference_composition
        ) {
            return active_solvation_factor_observation(
                active_model,
                active_parameters,
                temperature,
                pressure,
                reference_composition,
                solvation_factor
            );
        },
        checkpoint
    );
    if (!reference.derivative.has_value() || reference.derivative->size() != 2
        || !reference.derivative_convergence_error.has_value()) {
        throw std::runtime_error(
            "active solvation-factor derivative is unavailable"
        );
    }
    auto reference_log_fugacity = reference.log_fugacity;
    double reference_molality = reference.molality;
    double reference_convergence_error = reference.convergence_error;
    if (solvation_factor == parameters_.solvation_factor[0]) {
        const auto value_reference = active_model.electrolyte_reference_sequence(
            temperature,
            pressure,
            std::vector<double>{1.0, 0.0, 0.0},
            std::array<std::size_t, 2>{1, 2},
            std::array<int, 2>{1, 1},
            "liquid",
            budget
        );
        reference_log_fugacity = value_reference.log_fugacity;
        reference_molality = value_reference.reference_molality;
        reference_convergence_error = value_reference.convergence_error;
    }
    std::vector<NativeAqueousMiacSolvationFactor> results(molalities.size());
    detail::cppad_parallel_for(molalities.size(), [&](std::size_t index) {
        checkpoint();
        const double molality = molalities[index];
        const auto composition = detail::formula_composition(
            std::vector<double>{1.0, 0.0, 0.0},
            active_parameters.molar_mass_kg_per_mol,
            std::array<std::size_t, 2>{1, 2},
            std::array<int, 2>{1, 1},
            molality
        );
        const auto finite = active_solvation_factor_observation(
            active_model,
            active_parameters,
            temperature,
            pressure,
            composition,
            solvation_factor
        );
        if (!finite.derivative.has_value() || finite.derivative->size() != 2) {
            throw std::runtime_error(
                "active solvation-factor derivative is unavailable"
            );
        }
        auto finite_log_fugacity = finite.log_fugacity;
        if (solvation_factor == parameters_.solvation_factor[0]) {
            const auto finite_state = active_model.evaluate_pressure(
                temperature, pressure, composition, "liquid"
            );
            if (!finite_state.log_fugacity_coefficient.has_value()) {
                throw std::runtime_error("finite fugacity is unavailable");
            }
            finite_log_fugacity = {
                (*finite_state.log_fugacity_coefficient)[1],
                (*finite_state.log_fugacity_coefficient)[2],
            };
        }
        const double conversion = std::log(
            1.0 + 2.0 * active_parameters.molar_mass_kg_per_mol[0] * molality
        );
        results[index] = {
            0.5 * (
                finite_log_fugacity[0] + finite_log_fugacity[1]
                - reference_log_fugacity[0] - reference_log_fugacity[1]
            ) - conversion,
            0.5 * (
                (*finite.derivative)[0] + (*finite.derivative)[1]
                - (*reference.derivative)[0] - (*reference.derivative)[1]
            ),
            reference_molality,
            reference_convergence_error,
            *reference.derivative_convergence_error,
        };
        const auto& result = results[index];
        if (!std::isfinite(
                result.log_mean_ionic_activity_coefficient_molality
            )
            || !std::isfinite(result.derivative)
            || !std::isfinite(result.reference_molality_mol_per_kg)
            || !std::isfinite(result.reference_convergence_error)
            || !std::isfinite(
                result.reference_derivative_convergence_error
            )) {
            throw std::runtime_error(
                "active solvation-factor result contains nonfinite values"
            );
        }
    });
    return results;
}

// EqID: packing_volume_bounds
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
    // EqID: compressibility_factor
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
        // EqID: residual_chemical_potential
        chemical_potential[component] = residual + compressibility - 1.0
            + derivatives.composition[component]
            - weighted_composition_derivative;
        if (log_fugacity.has_value()) {
            // EqID: log_fugacity_coefficient
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
        // EqID: residual_heat_capacity_cp
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
        electrolyte.bulk_relative_permittivity,
        residual,
        compressibility,
        // EqID: pressure
        compressibility * molar_density * detail::GAS_CONSTANT * temperature,
        std::move(chemical_potential),
        std::move(log_fugacity),
        // EqID: residual_entropy
        detail::GAS_CONSTANT
            * (-residual - temperature * derivatives.temperature),
        // EqID: residual_internal_energy
        -detail::GAS_CONSTANT * temperature * temperature
            * derivatives.temperature,
        // EqID: residual_enthalpy
        detail::GAS_CONSTANT * temperature
            * (-temperature * derivatives.temperature + compressibility - 1.0),
        // EqID: residual_heat_capacity_cv
        detail::GAS_CONSTANT * cv_over_r,
        residual_cp,
        molar_density,
        // EqID: pressure_density_derivative
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
    // EqID: density_upper_bound
    const double packing_coefficient = detail::packing_fraction_coefficient(
        temperature, mole_fractions, parameters_
    );
    if (!(packing_coefficient > 0.0) || !std::isfinite(packing_coefficient)) {
        throw std::invalid_argument("invalid pressure-state packing coefficient");
    }
    const auto root = detail::solve_density(
        pressure,
        // EqID: ideal_gas_density
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
