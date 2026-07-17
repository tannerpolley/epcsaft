#include "eos.hpp"

#include "association.hpp"
#include "born.hpp"
#include "electrolyte.hpp"
#include "helmholtz.hpp"

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
        molar_density,
        mole_fractions,
        parameters_,
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

    using ADScalar = CppAD::AD<double>;
    const std::size_t sites = site_fraction.size();
    std::vector<ADScalar> independent(1 + sites);
    independent[0] = molar_density;
    for (std::size_t site = 0; site < sites; ++site) {
        independent[1 + site] = site_fraction[site];
    }
    CppAD::Independent(independent);
    std::vector<ADScalar> taped_fraction(sites);
    for (std::size_t site = 0; site < sites; ++site) {
        taped_fraction[site] = independent[1 + site];
    }
    const auto taped = detail::residual_helmholtz(
        temperature,
        independent[0],
        mole_fractions,
        parameters_,
        sigma,
        diameter,
        taped_fraction,
        detail::electrolyte_contribution(
            temperature,
            independent[0],
            molar_density,
            mole_fractions,
            parameters_,
            sigma
        ).debye_huckel,
        ADScalar(born_ssm_ds)
    );
    const auto delta = detail::association_delta(
        temperature,
        independent[0],
        mole_fractions,
        parameters_,
        sigma,
        diameter
    );
    const auto mass_action = detail::association_residuals(
        independent[0],
        mole_fractions,
        parameters_,
        taped_fraction,
        delta
    );
    std::vector<ADScalar> dependent = {
        taped.hard_chain + taped.dispersion + taped.association
            + taped.debye_huckel + taped.born_ssm_ds,
    };
    dependent.insert(dependent.end(), mass_action.begin(), mass_action.end());
    CppAD::ADFun<double> function(independent, dependent);
    std::vector<double> state = {molar_density};
    state.insert(state.end(), site_fraction.begin(), site_fraction.end());
    const auto jacobian = function.Jacobian(state);
    const std::size_t width = 1 + sites;
    double density_derivative = jacobian[0];
    if (sites != 0) {
        std::vector<double> h_x(sites * sites);
        std::vector<double> negative_h_rho(sites);
        for (std::size_t row = 0; row < sites; ++row) {
            negative_h_rho[row] = -jacobian[(1 + row) * width];
            for (std::size_t column = 0; column < sites; ++column) {
                h_x[row * sites + column] =
                    jacobian[(1 + row) * width + 1 + column];
            }
        }
        const auto fraction_derivative =
            detail::solve_linear(std::move(h_x), std::move(negative_h_rho));
        for (std::size_t site = 0; site < sites; ++site) {
            density_derivative += jacobian[1 + site] * fraction_derivative[site];
        }
    }

    const double residual =
        values.hard_chain + values.dispersion + values.association
        + values.debye_huckel + values.born_ssm_ds;
    const double compressibility = 1.0 + molar_density * density_derivative;
    return {
        values.hard_chain,
        values.dispersion,
        values.association,
        values.debye_huckel,
        values.born_ssm_ds,
        residual,
        compressibility,
        compressibility * molar_density * detail::GAS_CONSTANT * temperature,
    };
}

}  // namespace epcsaft
