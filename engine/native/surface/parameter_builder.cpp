#include "surface/parameter_builder.hpp"

#include "core/trace.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

namespace epcsaft {
namespace {

using Index = std::vector<std::size_t>;

[[noreturn]] void reject(const std::string& message) {
    throw std::invalid_argument(message);
}

std::string where(const ParameterEntry& entry) {
    return entry.identity.empty() ? entry.field : entry.identity;
}

bool symmetric_field(const std::string& field) {
    return field == "k_ij" || field == "l_ij"
        || field == "k_ij_reciprocal_temperature_slope"
        || field == "k_ij_reference_temperature"
        || field == "association_energy_over_k"
        || field == "association_volume";
}

Index canonical_index(const ParameterEntry& entry) {
    Index result = entry.index;
    if (symmetric_field(entry.field) && result.size() == 2 && result[0] > result[1]) {
        std::swap(result[0], result[1]);
    }
    return result;
}

double unit_factor(const ParameterEntry& entry, const std::string& expected) {
    if (expected == "molar mass") {
        // The packet vocabulary is deliberately closed and source-neutral.
        const std::array<std::pair<const char*, double>, 3> units = {
            std::pair{"kilogram / mole", 1.0},
            std::pair{"kg / mol", 1.0},
            std::pair{"gram / mole", 1.0e-3},
        };
        for (const auto& [unit, factor] : units) {
            if (entry.unit == unit) return factor;
        }
    } else if (expected == "dimensionless" && entry.unit == "dimensionless") {
        return 1.0;
    } else if ((expected == "kelvin" && entry.unit == "kelvin") || (expected == "length" && entry.unit == "angstrom")) {
        return 1.0;
    } else if (expected == "charge" && entry.unit == "elementary-charge") {
        return 1.0;
    } else if (expected == "reciprocal temperature" && entry.unit == "1 / kelvin") {
        return 1.0;
    } else if (expected == "molar volume" && entry.unit == "liter / mole") {
        return 1.0e-3;
    } else if ((expected == "debye" || expected == "debye ** 2" || expected == "angstrom ** 3") && entry.unit == expected) {
        return 1.0;
    }
    std::ostringstream message;
    message << where(entry) << ": expected " << expected;
    reject("unit " + entry.unit + " for " + message.str());
}

double unit_factor(const ParameterEntry& entry, int power) {
    const std::string expected = power == 1
        ? "1 / kelvin"
        : "1 / kelvin ** " + std::to_string(power);
    if (entry.unit != expected) {
        reject("unit " + entry.unit + " for " + where(entry) + ": expected " + expected);
    }
    return 1.0;
}

void require_index(const ParameterEntry& entry, std::size_t size, const char* kind) {
    if (entry.index.size() != size) {
        reject("invalid index for " + where(entry) + ": expected " + kind);
    }
}

std::size_t component_index(const ParameterEntry& entry, std::size_t components) {
    require_index(entry, 1, "component index");
    if (entry.index[0] >= components) reject("component index out of range for " + where(entry));
    return entry.index[0];
}

std::pair<std::size_t, std::size_t> pair_index(
    const ParameterEntry& entry,
    std::size_t components
) {
    require_index(entry, 2, "pair index");
    if (entry.index[0] >= components || entry.index[1] >= components
        || entry.index[0] == entry.index[1]) {
        reject("invalid component pair for " + where(entry));
    }
    return {entry.index[0], entry.index[1]};
}

std::pair<std::size_t, std::size_t> site_index(
    const ParameterEntry& entry,
    std::size_t sites
) {
    require_index(entry, 2, "site index");
    if (entry.index[0] >= sites || entry.index[1] >= sites) {
        reject("site index out of range for " + where(entry));
    }
    return {entry.index[0], entry.index[1]};
}

struct Rule final {
    std::string identity;
    std::array<std::size_t, 6> index{};
};

}  // namespace

Parameters build_parameters(const ParameterDocument& document) {
    core::Scope scope("build_parameters");
    const std::size_t components = document.component_ids.size();
    const std::size_t sites = document.site_component.size();
    if (document.site_multiplicity.size() != sites) reject("site structural lists have different sizes");
    if (document.base_choice != "pc-saft") reject("unknown base choice " + document.base_choice);
    if (document.association_choice == "none") {
        if (sites != 0) reject("association none has unconsumed site topology");
    } else if (document.association_choice != "general-site") {
        reject("unknown association choice " + document.association_choice);
    }

    Parameters result;
    result.segment_count.assign(components, 0.0);
    result.sigma_constant.assign(components, 0.0);
    result.dispersion_energy_over_k.assign(components, 0.0);
    result.molar_mass_kg_per_mol.assign(components, 0.0);
    result.charge_number.assign(components, 0);
    result.packing_diameter_angstrom.assign(components, 0.0);
    result.debye_huckel_diameter_angstrom.assign(components, 0.0);
    result.born_diameter_angstrom.assign(components, 0.0);
    result.solvation_factor.assign(components, 0.0);
    result.relative_permittivity.assign(components, 0.0);
    result.permittivity_reference_temperature_k.assign(components, 0.0);
    result.ion_suppression_coefficient.assign(components, 0.0);
    result.dielectric_volume_m3_per_mol.assign(components, 0.0);
    result.dielectric_temperature_k.assign(components, 0.0);
    result.dipole_scaling_debye2.assign(components, 0.0);
    result.polarizability_scaling_angstrom3.assign(components, 0.0);
    result.correlation_integral_parameter.assign(components, 0.0);
    result.k_ij.assign(components * components, 0.0);
    result.k_ij_reciprocal_temperature_slope.assign(components * components, 0.0);
    result.k_ij_reference_temperature_k.assign(components * components, 0.0);
    result.l_ij.assign(components * components, 0.0);
    result.site_component = document.site_component;
    result.site_multiplicity = document.site_multiplicity;
    result.association_active.assign(sites * sites, 0);
    result.association_energy_over_k.assign(sites * sites, 0.0);
    result.association_volume.assign(sites * sites, 0.0);

    if (document.electrolyte_choice == "fully-dissociated-debye-huckel") {
        result.electrolyte_form = ElectrolyteForm::DebyeHuckel;
    } else if (document.electrolyte_choice == "born") {
        result.electrolyte_form = ElectrolyteForm::DebyeHuckelBorn;
    } else if (document.electrolyte_choice != "none") {
        reject("unknown electrolyte choice " + document.electrolyte_choice);
    }
    if (document.electrolyte_choice == "none") {
        if (!document.ionic_dispersion_choice.empty()) {
            reject("unconsumed model/electrolyte/ionic_dispersion");
        }
    } else if (document.ionic_dispersion_choice.empty()) {
        reject("missing model/electrolyte/ionic_dispersion");
    } else if (document.ionic_dispersion_choice == "all-ion-pairs") {
        result.ionic_dispersion = IonicDispersion::AllIonPairs;
    } else if (document.ionic_dispersion_choice == "exclude-same-sign-ion-pairs") {
        result.ionic_dispersion = IonicDispersion::ExcludeSameSign;
    } else {
        reject("unknown ionic_dispersion choice " + document.ionic_dispersion_choice);
    }

    if (document.permittivity_choice == "solvent-only"
        || document.permittivity_choice == "ion-fraction-suppression"
        || document.permittivity_choice == "ion-specific-suppression") {
        result.permittivity_form = PermittivityForm::IonSuppression;
    } else if (document.permittivity_choice == "component-permittivity-mixing") {
        result.permittivity_form = PermittivityForm::ComponentMixing;
    } else if (document.permittivity_choice == "solvent-density-temperature") {
        result.permittivity_form = PermittivityForm::DensityTemperature;
    } else if (document.permittivity_choice == "dipolar-perturbation") {
        result.permittivity_form = PermittivityForm::DipolarPerturbation;
    } else if (document.permittivity_choice != "none") {
        reject("unknown permittivity choice " + document.permittivity_choice);
    }

    // The polar family is optional (ADR 0011): absent or none adds no term.
    const bool polar = document.polar_choice == "dipole-dipole";
    if (!polar && !document.polar_choice.empty() && document.polar_choice != "none") reject("unknown polar choice " + document.polar_choice);
    if (polar) {
        result.dipole_moment_debye.assign(components, 0.0);
        result.dipole_count.assign(components, 0.0);
    }

    std::set<std::tuple<std::string, Index>> seen;
    std::set<std::pair<std::string, std::size_t>> component_fields;
    std::set<std::string> scalar_fields;
    double ion_fraction_suppression = 0.0;
    std::vector<Rule> rules;
    std::map<std::size_t, std::size_t> sigma_term_components;
    std::map<std::size_t, std::size_t> permittivity_term_components;
    // Density-dependent forms consume only their own records (full EOS 6.2, 6.4).
    const bool density_temperature = document.permittivity_choice == "solvent-density-temperature";
    const bool dipolar = document.permittivity_choice == "dipolar-perturbation";
    for (const auto& entry : document.entries) {
        const bool dipolar_field = entry.field == "dipole_scaling" || entry.field == "polarizability_scaling"
            || entry.field == "correlation_integral_parameter";
        const bool density_field = entry.field == "dielectric_volume" || entry.field == "dielectric_temperature";
        const bool dipole_field = entry.field == "dipole_moment" || entry.field == "dipole_count";
        if ((entry.field == "ion_suppression_coefficient"
             && document.permittivity_choice != "ion-fraction-suppression"
             && document.permittivity_choice != "ion-specific-suppression")
            || (dipolar_field && !dipolar) || (density_field && !density_temperature) || (dipole_field && !polar)
            || ((entry.field == "relative_permittivity" || entry.field == "permittivity_term_coefficient"
                 || entry.field == "permittivity_reference_temperature") && (dipolar || density_temperature))
            || ((entry.field == "ionic_region_relative_permittivity"
                 || entry.field == "born_c_shell"
                 || entry.field == "born_c_dielectric"
                 || entry.field == "born_diameter"
                 || entry.field == "solvation_factor")
                && document.electrolyte_choice != "born")) {
            reject("unconsumed " + where(entry));
        }
        const auto key = std::make_tuple(entry.field, canonical_index(entry));
        if (!seen.emplace(key).second) reject("duplicate " + where(entry));
        double factor = 1.0;
        if (entry.field == "molar_mass") {
            factor = unit_factor(entry, "molar mass");
            const auto i = component_index(entry, components);
            component_fields.emplace(entry.field, i);
            result.molar_mass_kg_per_mol[i] = entry.magnitude * factor;
        } else if (entry.field == "charge_number") {
            factor = unit_factor(entry, "charge");
            const auto i = component_index(entry, components);
            component_fields.emplace(entry.field, i);
            if (!std::isfinite(entry.magnitude) || std::trunc(entry.magnitude) != entry.magnitude) {
                reject("charge " + where(entry) + " is not integral");
            }
            result.charge_number[i] = static_cast<int>(entry.magnitude * factor);
        } else if (entry.field == "segment_count" || entry.field == "sigma_constant"
                   || entry.field == "dispersion_energy_over_k" || entry.field == "solvation_factor"
                   || entry.field == "relative_permittivity" || entry.field == "packing_diameter"
                   || entry.field == "debye_huckel_diameter" || entry.field == "born_diameter"
                   || (entry.field == "ion_suppression_coefficient" && !entry.index.empty())) {
            const auto i = component_index(entry, components);
            component_fields.emplace(entry.field, i);
            const std::string dimension = entry.field == "segment_count"
                || entry.field == "solvation_factor" || entry.field == "relative_permittivity"
                || entry.field == "ion_suppression_coefficient" ? "dimensionless"
                : (entry.field == "dispersion_energy_over_k" ? "kelvin" : "length");
            factor = unit_factor(entry, dimension);
            if (entry.field == "segment_count") result.segment_count[i] = entry.magnitude * factor;
            else if (entry.field == "sigma_constant") result.sigma_constant[i] = entry.magnitude * factor;
            else if (entry.field == "dispersion_energy_over_k") result.dispersion_energy_over_k[i] = entry.magnitude * factor;
            else if (entry.field == "solvation_factor") result.solvation_factor[i] = entry.magnitude * factor;
            else if (entry.field == "relative_permittivity") result.relative_permittivity[i] = entry.magnitude * factor;
            else if (entry.field == "packing_diameter") result.packing_diameter_angstrom[i] = entry.magnitude * factor;
            else if (entry.field == "debye_huckel_diameter") result.debye_huckel_diameter_angstrom[i] = entry.magnitude * factor;
            else if (entry.field == "born_diameter") result.born_diameter_angstrom[i] = entry.magnitude * factor;
            else result.ion_suppression_coefficient[i] = entry.magnitude * factor;
        } else if (density_field || dipolar_field) {
            const auto i = component_index(entry, components);
            component_fields.emplace(entry.field, i);
            factor = unit_factor(entry, entry.field == "dielectric_volume" ? "molar volume"
                : entry.field == "dielectric_temperature" ? "kelvin"
                : entry.field == "dipole_scaling" ? "debye ** 2"
                : entry.field == "polarizability_scaling" ? "angstrom ** 3" : "dimensionless");
            const double value = entry.magnitude * factor;
            if (entry.field == "dielectric_volume") result.dielectric_volume_m3_per_mol[i] = value;
            else if (entry.field == "dielectric_temperature") result.dielectric_temperature_k[i] = value;
            else if (entry.field == "dipole_scaling") result.dipole_scaling_debye2[i] = value;
            else if (entry.field == "polarizability_scaling") result.polarizability_scaling_angstrom3[i] = value;
            else result.correlation_integral_parameter[i] = value;
        } else if (dipole_field) {
            const auto i = component_index(entry, components);
            component_fields.emplace(entry.field, i);
            factor = unit_factor(entry, entry.field == "dipole_moment" ? "debye" : "dimensionless");
            if (!std::isfinite(entry.magnitude) || entry.magnitude < 0.0) reject("negative or nonfinite " + where(entry));
            (entry.field == "dipole_moment" ? result.dipole_moment_debye : result.dipole_count)[i] = entry.magnitude * factor;
        } else if (entry.field == "sigma_term_amplitude" || entry.field == "sigma_term_exponent") {
            require_index(entry, 2, "sigma term index");
            if (entry.index[0] >= sigma_term_components.size() + 1 || entry.index[1] >= components) {
                reject("sigma term index out of range for " + where(entry));
            }
            sigma_term_components[entry.index[0]] = entry.index[1];  // units are converted in the term pass below
        } else if (entry.field == "permittivity_term_coefficient") {
            require_index(entry, 3, "permittivity term index");
            if (entry.index[0] >= permittivity_term_components.size() + 1 || entry.index[1] >= components) {
                reject("permittivity term index out of range for " + where(entry));
            }
            permittivity_term_components[entry.index[0]] = entry.index[1];
        } else if (entry.field == "permittivity_reference_temperature") {
            const auto i = component_index(entry, components);
            component_fields.emplace(entry.field, i);
            factor = unit_factor(entry, "kelvin");
            result.permittivity_reference_temperature_k[i] = entry.magnitude * factor;
        } else if (entry.field == "k_ij" || entry.field == "l_ij"
                   || entry.field == "k_ij_reciprocal_temperature_slope"
                   || entry.field == "k_ij_reference_temperature") {
            const auto [i, j] = pair_index(entry, components);
            const auto offset = i * components + j;
            const std::string dimension = entry.field == "k_ij" || entry.field == "l_ij"
                ? "dimensionless" : "kelvin";
            factor = unit_factor(entry, dimension);
            if (entry.field == "k_ij") {
                result.k_ij[offset] = entry.magnitude * factor;
                result.k_ij[j * components + i] = entry.magnitude * factor;
            } else if (entry.field == "l_ij") {
                result.l_ij[offset] = entry.magnitude * factor;
                result.l_ij[j * components + i] = entry.magnitude * factor;
            } else if (entry.field == "k_ij_reciprocal_temperature_slope") {
                result.k_ij_reciprocal_temperature_slope[offset] = entry.magnitude * factor;
                result.k_ij_reciprocal_temperature_slope[j * components + i] = entry.magnitude * factor;
            } else {
                result.k_ij_reference_temperature_k[offset] = entry.magnitude * factor;
                result.k_ij_reference_temperature_k[j * components + i] = entry.magnitude * factor;
            }
        } else if (entry.field == "association_energy_over_k" || entry.field == "association_volume") {
            const auto [s, t] = site_index(entry, sites);
            factor = unit_factor(entry, entry.field == "association_volume" ? "dimensionless" : "kelvin");
            result.association_active[s * sites + t] = 1;
            result.association_active[t * sites + s] = 1;
            if (entry.field == "association_energy_over_k") {
                result.association_energy_over_k[s * sites + t] = entry.magnitude * factor;
                result.association_energy_over_k[t * sites + s] = entry.magnitude * factor;
            } else {
                result.association_volume[s * sites + t] = entry.magnitude * factor;
                result.association_volume[t * sites + s] = entry.magnitude * factor;
            }
        } else if (entry.field == "association_combining_rule") {
            if (entry.identity != "arithmetic-energy-geometric-volume" && entry.identity != "wolbach-sandler") reject("unknown association rule " + entry.identity);
            require_index(entry, 6, "combining rule index");
            for (const auto index : entry.index) if (index >= sites) reject("site index out of range for " + where(entry));
            rules.push_back({entry.identity, {entry.index[0], entry.index[1], entry.index[2], entry.index[3], entry.index[4], entry.index[5]}});
        } else if (entry.field == "ionic_region_relative_permittivity"
                   || (entry.field == "ion_suppression_coefficient" && entry.index.empty())) {
            require_index(entry, 0, "scalar index");
            scalar_fields.insert(entry.field);
            factor = unit_factor(entry, "dimensionless");
            (entry.field == "ionic_region_relative_permittivity" ? result.ionic_region_relative_permittivity : ion_fraction_suppression)
                = entry.magnitude * factor;
        } else if (entry.field == "born_c_shell" || entry.field == "born_c_dielectric") {
            require_index(entry, 0, "scalar index");
            scalar_fields.insert(entry.field);
            factor = unit_factor(entry, "dimensionless");
            if (entry.field == "born_c_shell") result.born_c_shell = entry.magnitude * factor;
            else result.born_c_dielectric = entry.magnitude * factor;
        } else reject("unknown parameter field " + entry.field + " for " + where(entry));
    }

    const std::size_t sigma_terms = sigma_term_components.empty() ? 0 : sigma_term_components.rbegin()->first + 1;
    result.sigma_term_component.assign(sigma_terms, 0);
    result.sigma_term_amplitude.assign(sigma_terms, 0.0);
    result.sigma_term_exponent.assign(sigma_terms, 0.0);
    const std::size_t permittivity_terms = permittivity_term_components.empty() ? 0 : permittivity_term_components.rbegin()->first + 1;
    result.permittivity_term_component.assign(permittivity_terms, 0);
    result.permittivity_term_coefficient.assign(permittivity_terms, 0.0);
    result.permittivity_term_power.assign(permittivity_terms, 0);
    std::set<std::size_t> sigma_amplitude, sigma_exponent, permittivity_coefficients;
    for (const auto& entry : document.entries) {
        if (entry.field == "sigma_term_amplitude" || entry.field == "sigma_term_exponent") {
            const auto term = entry.index[0];
            result.sigma_term_component[term] = static_cast<int>(entry.index[1]);
            const double value = entry.field == "sigma_term_amplitude"
                ? entry.magnitude * unit_factor(entry, "length")
                : entry.magnitude * unit_factor(entry, "reciprocal temperature");
            if (entry.field == "sigma_term_amplitude") {
                result.sigma_term_amplitude[term] = value; sigma_amplitude.insert(term);
            } else { result.sigma_term_exponent[term] = value; sigma_exponent.insert(term); }
        } else if (entry.field == "permittivity_term_coefficient") {
            const auto term = entry.index[0];
            result.permittivity_term_component[term] = static_cast<int>(entry.index[1]);
            result.permittivity_term_power[term] = static_cast<int>(entry.index[2]);
            result.permittivity_term_coefficient[term] = entry.magnitude * unit_factor(entry, static_cast<int>(entry.index[2]));
            permittivity_coefficients.insert(term);
        }
    }
    if (sigma_amplitude != sigma_exponent || sigma_amplitude.size() != sigma_terms) reject("incomplete sigma term");
    if (permittivity_coefficients.size() != permittivity_terms) reject("incomplete permittivity term");

    for (std::size_t i = 0; i < components; ++i) {
        if (component_fields.find({"charge_number", i}) == component_fields.end()) reject("missing component/" + document.component_ids[i] + "/charge_number");
        if (result.charge_number[i] == 0 && (component_fields.find({"molar_mass", i}) == component_fields.end()
                                             || result.molar_mass_kg_per_mol[i] <= 0.0)) {
            reject("missing component/" + document.component_ids[i] + "/molar_mass");
        }
        if (component_fields.find({"segment_count", i}) == component_fields.end()) reject("missing component/" + document.component_ids[i] + "/segment_count");
        if (component_fields.find({"sigma_constant", i}) == component_fields.end()) reject("missing component/" + document.component_ids[i] + "/segment_diameter");
        if (component_fields.find({"dispersion_energy_over_k", i}) == component_fields.end()) reject("missing component/" + document.component_ids[i] + "/dispersion_energy_over_k");
        if (result.electrolyte_form.is(ElectrolyteForm::DebyeHuckelBorn) && component_fields.find({"solvation_factor", i}) == component_fields.end()) reject("missing component/" + document.component_ids[i] + "/solvation_factor");
        if (result.charge_number[i] != 0 && component_fields.find({"packing_diameter", i}) == component_fields.end()) reject("missing component/" + document.component_ids[i] + "/packing_diameter");
        if (result.charge_number[i] != 0 && component_fields.find({"debye_huckel_diameter", i}) == component_fields.end()) reject("missing component/" + document.component_ids[i] + "/debye_huckel_diameter");
        if (result.electrolyte_form.is(ElectrolyteForm::DebyeHuckelBorn) && result.charge_number[i] != 0 && component_fields.find({"born_diameter", i}) == component_fields.end()) reject("missing component/" + document.component_ids[i] + "/born_diameter");
        // Schreckenberg records belong to every neutral solvent only; Rueben records to every component.
        for (const char* field : {"dielectric_volume", "dielectric_temperature"}) {
            const bool present = component_fields.find({field, i}) != component_fields.end();
            if (density_temperature && present != (result.charge_number[i] == 0)) reject((present ? "unconsumed component/" : "missing component/") + document.component_ids[i] + "/" + field);
        }
        for (const char* field : {"dipole_scaling", "polarizability_scaling", "correlation_integral_parameter"}) {
            if (dipolar && component_fields.find({field, i}) == component_fields.end()) reject("missing component/" + document.component_ids[i] + "/" + field);
        }
        for (const char* field : {"dipole_moment", "dipole_count"}) {
            if (polar && component_fields.find({field, i}) == component_fields.end()) reject("missing component/" + document.component_ids[i] + "/" + field);
        }
    }
    for (std::size_t i = 0; i < components; ++i) for (std::size_t j = i + 1; j < components; ++j) {
        const bool same_sign = result.charge_number[i] * result.charge_number[j] > 0;
        if (same_sign && result.ionic_dispersion.is(IonicDispersion::ExcludeSameSign)) {
            for (const auto& entry : document.entries) {
                if (entry.index.size() != 2
                    || !((entry.index[0] == i && entry.index[1] == j)
                         || (entry.index[0] == j && entry.index[1] == i))) {
                    continue;
                }
                if (entry.field == "k_ij"
                    || entry.field == "k_ij_reciprocal_temperature_slope"
                    || entry.field == "k_ij_reference_temperature") {
                    const std::string identity = where(entry);
                    reject("unconsumed " + identity + " under exclude-same-sign-ion-pairs");
                }
            }
            continue;
        }
        if (result.k_ij[i * components + j] == 0.0 && result.k_ij[j * components + i] == 0.0) {
            bool present = false;
            for (const auto& entry : document.entries) if (entry.field == "k_ij" && entry.index.size() == 2 && ((entry.index[0] == i && entry.index[1] == j) || (entry.index[0] == j && entry.index[1] == i))) present = true;
            if (!present) {
                const auto& first = std::min(document.component_ids[i], document.component_ids[j]);
                const auto& second = std::max(document.component_ids[i], document.component_ids[j]);
                reject("missing pair/" + first + "/" + second + "/k_ij");
            }
        }
        if (result.k_ij_reciprocal_temperature_slope[i * components + j] != 0.0 && result.k_ij_reference_temperature_k[i * components + j] == 0.0) reject("missing pair reference temperature/" + document.component_ids[i] + "/" + document.component_ids[j]);
    }
    if (result.electrolyte_form.is(ElectrolyteForm::DebyeHuckelBorn)
        && scalar_fields.find("ionic_region_relative_permittivity") == scalar_fields.end()) reject("missing ionic_region_relative_permittivity");
    if (result.electrolyte_form.is(ElectrolyteForm::DebyeHuckelBorn)
        && scalar_fields.find("born_c_shell") == scalar_fields.end()) reject("missing born_c_shell");
    if (result.electrolyte_form.is(ElectrolyteForm::DebyeHuckelBorn)
        && scalar_fields.find("born_c_dielectric") == scalar_fields.end()) reject("missing born_c_dielectric");
    if (document.permittivity_choice == "ion-fraction-suppression") {
        if (!scalar_fields.count("ion_suppression_coefficient")) reject("missing model/relative_permittivity/ion_fraction_suppression_coefficient");
        for (std::size_t i = 0; i < components; ++i) result.ion_suppression_coefficient[i] = result.charge_number[i] == 0 ? 0.0 : ion_fraction_suppression;
    } else if (document.permittivity_choice == "solvent-only" || document.permittivity_choice == "none") {
        std::fill(result.ion_suppression_coefficient.begin(), result.ion_suppression_coefficient.end(), 0.0);
    } else if (document.permittivity_choice == "ion-specific-suppression") {
        for (std::size_t i = 0; i < components; ++i) if (result.charge_number[i] != 0 && component_fields.find({"ion_suppression_coefficient", i}) == component_fields.end()) reject("missing component/" + document.component_ids[i] + "/ion_suppression_coefficient");
    }
    for (const auto& rule : rules) {
        const auto output = rule.index[0] * sites + rule.index[1];
        const auto self1 = rule.index[2] * sites + rule.index[3];
        const auto self2 = rule.index[4] * sites + rule.index[5];
        if (result.association_active[self1] == 0 || result.association_active[self2] == 0) reject("combining rule self edge inactive: " + rule.identity);
        if (result.association_active[output] != 0) reject("duplicate " + rule.identity);
        const double energy = 0.5 * (result.association_energy_over_k[self1] + result.association_energy_over_k[self2]);
        const double volume = std::sqrt(result.association_volume[self1] * result.association_volume[self2]);
        const auto reverse = rule.index[1] * sites + rule.index[0];
        result.association_active[output] = result.association_active[reverse] = rule.identity == "wolbach-sandler" ? 2 : 1;
        result.association_energy_over_k[output] = energy;
        result.association_energy_over_k[reverse] = energy;
        result.association_volume[output] = volume;
        result.association_volume[reverse] = volume;
    }
    // Selection fixes topology (full EOS §12): an uncharged selection has no
    // ionic term, so its packet-level electrolyte choices have no consumer.
    if (std::none_of(result.charge_number.begin(), result.charge_number.end(), [](int z) { return z != 0; })) {
        result.electrolyte_form = {};
        result.ionic_dispersion = {};
        result.permittivity_form = {};
    }
    return result;
}

}  // namespace epcsaft
