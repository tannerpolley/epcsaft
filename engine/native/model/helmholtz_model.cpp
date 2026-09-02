#include "model/helmholtz_model.hpp"

#include "physics/association.hpp"
#include "physics/born.hpp"
#include "physics/electrolyte.hpp"
#include "physics/helmholtz.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
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

void require(bool condition, const char* field) {
    if (!condition) {
        throw std::invalid_argument(std::string("invalid ") + field);
    }
}

}  // namespace

HelmholtzModel::HelmholtzModel(
    Parameters parameters,
    std::vector<detail::ActiveParameter> active
)
    : parameters_(std::move(parameters)), active_(std::move(active)) {
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
    require(positive_vector(parameters_.molar_mass_kg_per_mol), "molar_mass");
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
        if (parameters_.electrolyte_form == ElectrolyteForm::DebyeHuckelBorn) {
            require(parameters_.born_diameter_angstrom[component] > 0.0,
                    "born_diameter");
        }
    }
    if (parameters_.electrolyte_form == ElectrolyteForm::DebyeHuckelBorn) {
        require(parameters_.ionic_region_relative_permittivity > 0.0,
                "ionic_region_relative_permittivity");
    }
    require(std::isfinite(parameters_.ionic_region_relative_permittivity),
            "ionic_region_relative_permittivity");
    require(std::isfinite(parameters_.born_c_shell), "born_c_shell");
    require(std::isfinite(parameters_.born_c_dielectric), "born_c_dielectric");
    if (parameters_.permittivity_form != PermittivityForm::None) {
        require(std::any_of(
            parameters_.charge_number.begin(), parameters_.charge_number.end(),
            [&](int charge) { return charge == 0; }
        ), "dielectric solvent");
    }

    std::set<std::pair<ParameterFamily, std::size_t>> seen;
    for (const auto& item : active_) {
        require(!item.indices.empty(), "active_parameter.indices");
        for (const std::size_t index : item.indices) {
            std::size_t size = 0;
            switch (item.family) {
            case ParameterFamily::SegmentCount: size = components; break;
            case ParameterFamily::SigmaConstant: size = components; break;
            case ParameterFamily::DispersionEnergyOverK: size = components; break;
            case ParameterFamily::Kij: size = pairs; break;
            case ParameterFamily::AssociationEnergyOverK: size = site_pairs; break;
            case ParameterFamily::AssociationVolume: size = site_pairs; break;
            }
            require(index < size, "active_parameter.indices");
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
    const std::vector<double>& composition
) const {
    detail::PredicateSignature result;
    double charge_square_sum = 0.0;
    for (std::size_t component = 0; component < composition.size(); ++component) {
        if (!dielectric_solvent_component(component, parameters_)) {
            const double charge = parameters_.charge_number[component];
            charge_square_sum += composition[component] * charge * charge;
        }
    }
    result.push_back(charge_square_sum == 0.0 ? 1 : 0);
    if (charge_square_sum == 0.0) return result;
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
                < policy::ELECTROLYTE_SERIES_THRESHOLD.value ? 1 : 0
        );
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
    validate_point(temperature, volume, amounts);
    const std::size_t components = parameters_.component_count();
    const std::size_t total_coordinates = 2 + components + active_.size();
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
    const auto electro = electrolyte_contribution(
        temperature, density, temperature, density, composition, composition,
        parameters_
    );
    const auto born = born_contribution(
        temperature, composition, parameters_, electro.bulk_relative_permittivity
    );
    const auto direct = residual_helmholtz(
        temperature, density, composition, parameters_, sigma, diameter, sites,
        electro.debye_huckel, born
    );
    const double direct_value = total * (
        direct.hard_chain + direct.dispersion + direct.association
        + direct.debye_huckel + direct.born
    );
    std::vector<double> point = {temperature, volume};
    point.insert(point.end(), amounts.begin(), amounts.end());
    for (const auto& item : active_) {
        const std::vector<double>* values = nullptr;
        switch (item.family) {
        case ParameterFamily::SegmentCount: values = &parameters_.segment_count; break;
        case ParameterFamily::SigmaConstant: values = &parameters_.sigma_constant; break;
        case ParameterFamily::DispersionEnergyOverK:
            values = &parameters_.dispersion_energy_over_k; break;
        case ParameterFamily::Kij: values = &parameters_.k_ij; break;
        case ParameterFamily::AssociationEnergyOverK:
            values = &parameters_.association_energy_over_k; break;
        case ParameterFamily::AssociationVolume: values = &parameters_.association_volume; break;
        }
        point.push_back((*values)[item.indices.front()]);
    }
    const auto signature = predicate_signature(temperature, density, composition);
    if (tape_ == nullptr || tape_->predicates != signature) {
        tape_ = detail::HelmholtzTape::record(
            parameters_, active_, point, sites, signature
        );
        ++counters_.recordings;
    }
    tape_->function.new_dynamic(sites);
    std::vector<double> weight(1 + sites.size(), 0.0);
    weight.front() = 1.0;
    HelmholtzDerivatives result;
    result.site_fractions = sites;
    const auto value = tape_->function.Forward(0, point);
    ++counters_.forward;
    result.value = value.front();
    if (order == DerivativeOrder::Value) return result;
    if (order == DerivativeOrder::Gradient) {
        result.gradient = tape_->function.Reverse(1, weight);
        ++counters_.reverse;
        return result;
    }
    if (order == DerivativeOrder::PressureSlope) {
        result.gradient = tape_->function.Reverse(1, weight);
        ++counters_.reverse;
        std::vector<double> direction(total_coordinates, 0.0);
        direction[1] = 1.0;
        const auto first = tape_->function.Forward(1, direction);
        ++counters_.forward;
        const auto second = tape_->function.Reverse(2, weight);
        ++counters_.reverse;
        result.volume_second_derivative = second[2 * 1 + 1];
        return result;
    }

    ++counters_.hessian_calls;
    result.gradient.resize(total_coordinates);
    result.hessian.assign(total_coordinates * total_coordinates, 0.0);
    const std::size_t sites_count = sites.size();
    Eigen::MatrixXd association_hessian = Eigen::MatrixXd::Zero(
        static_cast<Eigen::Index>(sites_count), static_cast<Eigen::Index>(sites_count)
    );
    const auto delta = association_delta(
        temperature, density, composition, parameters_, sigma, diameter
    );
    const double number_density = density * NUMBER_DENSITY_PER_ANGSTROM3;
    for (std::size_t a = 0; a < sites_count; ++a) {
        const std::size_t component_a = static_cast<std::size_t>(
            parameters_.site_component[a]
        );
        const double mass_a = composition[component_a]
            * parameters_.site_multiplicity[a];
        for (std::size_t b = 0; b < sites_count; ++b) {
            const std::size_t component_b = static_cast<std::size_t>(
                parameters_.site_component[b]
            );
            const double mass_b = composition[component_b]
                * parameters_.site_multiplicity[b];
            association_hessian(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)) =
                -total * (
                    (a == b ? mass_a / (sites[a] * sites[a]) : 0.0)
                    + number_density * mass_a * mass_b * delta[a * sites_count + b]
                );
        }
    }
    Eigen::FullPivLU<Eigen::MatrixXd> factorization(association_hessian);
    if (sites_count != 0 && factorization.rank() != static_cast<Eigen::Index>(sites_count)) {
        throw std::runtime_error("association Hessian is singular");
    }
    Eigen::MatrixXd cross = Eigen::MatrixXd::Zero(
        static_cast<Eigen::Index>(sites_count), static_cast<Eigen::Index>(total_coordinates)
    );
    Eigen::MatrixXd base = Eigen::MatrixXd::Zero(
        static_cast<Eigen::Index>(total_coordinates), static_cast<Eigen::Index>(total_coordinates)
    );
    for (std::size_t column = 0; column < total_coordinates; ++column) {
        std::vector<double> direction(total_coordinates, 0.0);
        direction[column] = 1.0;
        const auto first = tape_->function.Forward(1, direction);
        ++counters_.forward;
        for (std::size_t site = 0; site < sites_count; ++site) {
            cross(static_cast<Eigen::Index>(site), static_cast<Eigen::Index>(column)) =
                first[1 + site];
        }
        const auto second = tape_->function.Reverse(2, weight);
        ++counters_.reverse;
        for (std::size_t row = 0; row < total_coordinates; ++row) {
            base(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) =
                second[2 * row + 1];
            if (column == 0) result.gradient[row] = second[2 * row];
        }
    }
    Eigen::MatrixXd correction = Eigen::MatrixXd::Zero(
        static_cast<Eigen::Index>(total_coordinates),
        static_cast<Eigen::Index>(total_coordinates)
    );
    if (sites_count != 0) {
        correction = cross.transpose() * factorization.solve(cross);
    }
    for (std::size_t row = 0; row < total_coordinates; ++row) {
        for (std::size_t column = 0; column < total_coordinates; ++column) {
            result.hessian[row * total_coordinates + column] =
                base(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column))
                - correction(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column));
        }
    }
    result.volume_second_derivative = result.hessian[1 * total_coordinates + 1];
    (void)direct_value;
    return result;
}

State HelmholtzModel::state(
    double temperature,
    double molar_density,
    const std::vector<double>& mole_fractions
) const {
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
    result.residual_internal_energy = -GAS_CONSTANT_J_PER_MOL_K * temperature * phi_t;
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
    const auto root = solve_density(
        *this, temperature, pressure, mole_fractions, branch, anchor_density
    );
    return state(temperature, root.molar_density, mole_fractions);
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
    std::lock_guard<std::mutex> lock(mutex_);
    return counters_;
}

void HelmholtzModel::reset_counters() const {
    std::lock_guard<std::mutex> lock(mutex_);
    counters_ = {};
}

}  // namespace epcsaft
