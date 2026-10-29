#include "surface/parameter_builder.hpp"

#include "core/trace.hpp"
#include "properties/helmholtz_model.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstring>
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
        || field == "association_volume"
        || field == "association_volume_correction" || field == "no_cross_association";
}

Index canonical_index(const ParameterEntry& entry) {
    Index result = entry.index;
    if (entry.field.rfind("cation_", 0) == 0 && result.size() == 4 && result[1] > result[2])
        std::swap(result[1], result[2]);
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
    } else if ((expected == "debye" || expected == "debye ** 2" || expected == "debye * angstrom" || expected == "angstrom ** 3"
                || expected == "kilogram / mole" || expected == "kilogram ** 2 / mole ** 2")
               && entry.unit == expected) {
        return 1.0;
    }
    std::ostringstream message;
    message << where(entry) << ": expected " << expected;
    reject("unit " + entry.unit + " for " + message.str());
}

// The unit kelvin ** exponent, in the closed vocabulary dimensionless, kelvin ** n and 1 / kelvin ** n.
double unit_factor(const ParameterEntry& entry, int exponent) {
    const int magnitude = std::abs(exponent);
    const std::string kelvin = magnitude == 1 ? "kelvin" : "kelvin ** " + std::to_string(magnitude);
    const std::string expected = exponent == 0 ? "dimensionless" : exponent > 0 ? kelvin : "1 / " + kelvin;
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

// Law targets in TemperatureLaw order; the volume correction (5) is consumed by a combining rule.
int law_target(const std::string& field) {
    const std::array<const char*, 5> targets = {"k_ij", "l_ij", "association_energy_over_k", "association_volume",
                                                "association_volume_correction"};
    for (std::size_t i = 0; i < targets.size(); ++i) if (field == targets[i]) return static_cast<int>(i) + 1;
    reject("unknown temperature law target " + field);
}

std::size_t canonical_slot(std::size_t a, std::size_t b, std::size_t size) {
    return std::min(a, b) * size + std::max(a, b);
}

}  // namespace

Parameters build_parameters(const ParameterDocument& document) {
    core::Scope scope("build_parameters");
    const std::size_t components = document.component_ids.size();
    const std::size_t sites = document.site_component.size();
    if (document.site_multiplicity.size() != sites) reject("site structural lists have different sizes");
    if (document.base_choice == "ideal-gas") {
        if (sites || !document.ionic_dispersion_choice.empty() || document.association_choice != "none" || document.electrolyte_choice != "none"
            || document.permittivity_choice != "none" || (!document.polar_choice.empty() && document.polar_choice != "none"))
            reject("ideal-gas has unconsumed residual families or sites");
        Parameters ideal;
        ideal.ideal_gas = true;
        ideal.molar_mass_kg_per_mol.assign(components, 0.0);
        ideal.charge_number.assign(components, 0);
        std::set<std::pair<std::string, Index>> seen;
        for (const auto& entry : document.entries) {
            if (!seen.emplace(entry.field, entry.index).second) reject("repeated " + where(entry));
            if (entry.field != "molar_mass" && entry.field != "charge_number") reject("ideal-gas unconsumed " + where(entry));
            const auto i = component_index(entry, components);
            if (entry.field == "charge_number") {
                if (entry.magnitude != 0.0 || unit_factor(entry, "dimensionless") != 1.0) reject("ideal-gas requires neutral species");
            } else ideal.molar_mass_kg_per_mol[i] = entry.magnitude * unit_factor(entry, "molar mass");
        }
        if (components == 0 || std::any_of(ideal.molar_mass_kg_per_mol.begin(), ideal.molar_mass_kg_per_mol.end(),
            [](double mass) { return !std::isfinite(mass) || mass <= 0.0; })) reject("ideal-gas requires positive molar masses");
        return ideal;
    }
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
    result.molecular_ion.assign(components, 0);
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
    } else if (document.permittivity_choice == "solvent-mole-weighting") {
        result.permittivity_form = PermittivityForm::SolventMoleWeighting;
        result.permittivity_solvent.assign(components, 0);
    } else if (document.permittivity_choice == "component-permittivity-mixing") {
        result.permittivity_form = PermittivityForm::ComponentMixing;
    } else if (document.permittivity_choice == "solvent-density-temperature") {
        result.permittivity_form = PermittivityForm::DensityTemperature;
    } else if (document.permittivity_choice == "dipolar-perturbation") {
        result.permittivity_form = PermittivityForm::DipolarPerturbation;
    } else if (document.permittivity_choice != "none") {
        reject("unknown permittivity choice " + document.permittivity_choice);
    }

    // The polar family is optional (ADRs 0011, 0014): absent or none adds no term.
    const bool polar = document.polar_choice == "dipole-dipole";
    const bool quadrupolar = document.polar_choice == "quadrupole-quadrupole";
    if (!polar && !quadrupolar && !document.polar_choice.empty() && document.polar_choice != "none") reject("unknown polar choice " + document.polar_choice);
    if (polar) {
        result.dipole_moment_debye.assign(components, 0.0);
        result.dipole_count.assign(components, 0.0);
    }
    if (quadrupolar) {
        result.quadrupole_moment_debye_angstrom.assign(components, 0.0);
        result.quadrupole_count.assign(components, 0.0);
    }

    std::set<std::tuple<std::string, Index>> seen;
    std::set<std::pair<std::string, std::size_t>> component_fields;
    std::set<std::string> scalar_fields;
    std::set<std::string> mixing_identities;
    double ion_fraction_suppression = 0.0;
    std::vector<Rule> rules;
    std::set<std::pair<std::size_t, std::size_t>> no_cross;
    std::map<std::size_t, std::size_t> sigma_term_components;
    std::map<std::size_t, std::size_t> permittivity_term_components;
    std::map<std::size_t, std::pair<std::string, double>> laws;  // quantity, reference (0: constant correction)
    std::set<std::pair<int, std::size_t>> law_slots;   // (target, canonical slot)
    std::map<std::size_t, std::size_t> corrections;   // canonical rule edge -> law
    // Density-dependent forms consume only their own records (full EOS 6.2, 6.4).
    const bool density_temperature = document.permittivity_choice == "solvent-density-temperature";
    const bool dipolar = document.permittivity_choice == "dipolar-perturbation";
    for (const auto& entry : document.entries) {
        const bool dipolar_field = entry.field == "dipole_scaling" || entry.field == "polarizability_scaling"
            || entry.field == "correlation_integral_parameter";
        const bool density_field = entry.field == "dielectric_volume" || entry.field == "dielectric_temperature";
        const bool dipole_field = entry.field == "dipole_moment" || entry.field == "dipole_count";
        const bool quadrupole_field = entry.field == "quadrupole_moment" || entry.field == "quadrupole_count";
        if ((entry.field == "ion_suppression_coefficient"
             && document.permittivity_choice != "ion-fraction-suppression"
             && document.permittivity_choice != "ion-specific-suppression")
            || (dipolar_field && !dipolar) || (density_field && !density_temperature) || (dipole_field && !polar) || (quadrupole_field && !quadrupolar)
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
        // A law [law, a, b] or volume correction owns its quantity's (a, b) value: one of value or law.
        const std::string suffix = "_temperature_law";
        const bool law = entry.field.size() > suffix.size() && entry.field.substr(entry.field.size() - suffix.size()) == suffix;
        const bool correction = law || entry.field == "association_volume_correction";
        ParameterEntry owned = entry;
        if (correction) {
            require_index(entry, 3, "temperature law index");
            owned.field = law ? entry.field.substr(0, entry.field.size() - suffix.size()) : entry.field;
            owned.index = {entry.index[1], entry.index[2]};
        }
        const auto key = std::make_tuple(owned.field, canonical_index(owned));
        if (!seen.emplace(key).second) reject("duplicate " + where(entry));
        double factor = 1.0;
        if (entry.field == "cation_pair_mixing" || entry.field == "cation_triple_mixing") {
            require_index(entry, 4, "solvent/cation/cation/anion indices");
            std::array<int, 4> topology{};
            const auto indices = canonical_index(entry);
            for (std::size_t r = 0; r < 4; ++r) {
                if (indices[r] >= components) reject("cation mixing index for " + where(entry));
                topology[r] = static_cast<int>(indices[r]);
            }
            if (entry.identity.empty() || !mixing_identities.insert(entry.identity).second) reject("cation mixing identity");
            if (result.cation_mixing.empty()) result.cation_mixing.push_back(topology);
            if (result.cation_mixing.front() != topology) reject("conflicting cation mixing topology");
            const bool pair = entry.field == "cation_pair_mixing";
            auto& values = pair ? result.cation_pair_mixing : result.cation_triple_mixing;
            values.push_back(entry.magnitude * unit_factor(entry, pair ? "kilogram / mole" : "kilogram ** 2 / mole ** 2"));
        } else if (entry.field == "molar_mass") {
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
        } else if (entry.field == "permittivity_solvent" || entry.field == "permittivity_normalization") {
            if (!result.permittivity_form.is(PermittivityForm::SolventMoleWeighting)) reject("unconsumed " + where(entry));
            const double value = entry.magnitude * unit_factor(entry, "dimensionless");
            if (entry.field == "permittivity_solvent") {
                const auto i = component_index(entry, components);
                if (value != 1.0) reject("invalid permittivity solvent membership");
                result.permittivity_solvent[i] = 1;
            } else {
                require_index(entry, 0, "permittivity normalization");
                if (value != 0.0 && value != 1.0) reject("invalid permittivity normalization");
                scalar_fields.insert(entry.field);
                result.permittivity_solvent_normalized = value == 1.0;
            }
        } else if (entry.field == "charged_molecule_ion") {
            const auto i = component_index(entry, components);
            if (entry.magnitude * unit_factor(entry, "dimensionless") != 1.0) reject("charged_molecule_ion is not 1: " + where(entry));
            result.molecular_ion[i] = 1;
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
        } else if (dipole_field || quadrupole_field) {
            const auto i = component_index(entry, components);
            component_fields.emplace(entry.field, i);
            factor = unit_factor(entry, entry.field == "dipole_moment" ? "debye"
                : entry.field == "quadrupole_moment" ? "debye * angstrom" : "dimensionless");
            if (!std::isfinite(entry.magnitude) || entry.magnitude < 0.0) reject("negative or nonfinite " + where(entry));
            (entry.field == "dipole_moment" ? result.dipole_moment_debye : entry.field == "dipole_count" ? result.dipole_count
                : entry.field == "quadrupole_moment" ? result.quadrupole_moment_debye_angstrom : result.quadrupole_count)[i] = entry.magnitude * factor;
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
        } else if (correction) {
            const int target = law_target(owned.field);
            const bool site = target > 2;
            const auto [a, b] = site ? site_index(owned, sites) : pair_index(owned, components);
            const std::pair declared{owned.field, law ? entry.magnitude * unit_factor(entry, "kelvin") : 0.0};
            if (law && !(std::isfinite(declared.second) && declared.second > 0.0)) reject("invalid reference temperature for " + where(entry));
            if (laws.emplace(entry.index[0], declared).first->second != declared) reject("conflicting temperature law " + where(entry));
            const auto slot = canonical_slot(a, b, site ? sites : components);
            if (target == 5) {
                corrections[slot] = entry.index[0];
            } else {
                law_slots.insert({target, slot});
                result.temperature_law_slot.push_back({target, static_cast<int>(slot), static_cast<int>(entry.index[0])});
                if (site) result.association_active[a * sites + b] = result.association_active[b * sites + a] = 1;
            }
        } else if (entry.field == "temperature_law_coefficient" || entry.field == "temperature_law_reciprocal_slope") {
            // Read in the term pass below, after every law is declared.
        } else if (entry.field == "no_cross_association") {
            // A declaration marker, not an association quantity: every unlike site edge stays inactive.
            const auto [i, j] = pair_index(entry, components);
            if (unit_factor(entry, "dimensionless") * entry.magnitude != 0.0) reject("nonzero declaration " + where(entry));
            no_cross.emplace(std::min(i, j), std::max(i, j));
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
    std::set<std::size_t> sigma_amplitude, sigma_exponent, permittivity_coefficients, intercepts;
    std::set<std::string> coefficient_identities;
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
            result.permittivity_term_coefficient[term] = entry.magnitude * unit_factor(entry, -static_cast<int>(entry.index[2]));
            permittivity_coefficients.insert(term);
        } else if (entry.field == "temperature_law_coefficient" || entry.field == "temperature_law_reciprocal_slope"
                   || (entry.field == "association_volume_correction" && !intercepts.count(entry.index[0]))) {
            // Coefficient a_p has unit K^(e-p) and s unit K^(e+1), e = 1 for an energy law and 0 otherwise.
            const bool slope = entry.field == "temperature_law_reciprocal_slope", constant = entry.field == "association_volume_correction";
            require_index(entry, slope ? 1 : (constant ? 3 : 2), "temperature law term index");
            const auto found = laws.find(entry.index[0]);
            if (found == laws.end()) reject("temperature law term without its law: " + where(entry));
            const int power = slope ? -1 : (constant ? 0 : static_cast<int>(entry.index[1]));
            if (power > 3 || (found->second.second == 0.0) != constant) reject("invalid temperature law term " + where(entry));
            const int energy = found->second.first == "association_energy_over_k" ? 1 : 0;
            const double value = entry.magnitude * unit_factor(entry, energy + (slope ? 1 : -power));
            if (!std::isfinite(value)) reject("nonfinite " + where(entry));
            if (!coefficient_identities.insert(entry.identity).second) reject("repeated temperature law coefficient " + where(entry));
            if (power == 0) intercepts.insert(entry.index[0]);
            result.temperature_law_term_law.push_back(static_cast<int>(entry.index[0]));
            result.temperature_law_term_power.push_back(power);
            result.temperature_law_coefficient_identity.push_back(entry.identity);
            result.temperature_law_coefficient.push_back(value);
        }
    }
    for (const auto& [number, law] : laws) {
        if (number != result.temperature_law_reference_k.size()) reject("temperature laws are not numbered consecutively");
        if (!intercepts.count(number)) reject("temperature law without its power-0 coefficient: " + law.first);
        result.temperature_law_reference_k.push_back(law.second);
    }
    if (sigma_amplitude != sigma_exponent || sigma_amplitude.size() != sigma_terms) reject("incomplete sigma term");
    if (permittivity_coefficients.size() != permittivity_terms) reject("incomplete permittivity term");

    if (result.permittivity_form.is(PermittivityForm::SolventMoleWeighting)) {
        if (!scalar_fields.count("permittivity_normalization")) reject("missing permittivity normalization");
        for (std::size_t i = 0; i < components; ++i) {
            const bool law = component_fields.count({"relative_permittivity", i})
                || std::find(result.permittivity_term_component.begin(), result.permittivity_term_component.end(), i) != result.permittivity_term_component.end();
            if (law != bool(result.permittivity_solvent[i])) reject("missing or unconsumed solvent dielectric law");
        }
    }
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
        if (result.molecular_ion[i] && result.charge_number[i] == 0) reject("unconsumed component/" + document.component_ids[i] + "/charged_molecule_ion");
        const bool packing = component_fields.find({"packing_diameter", i}) != component_fields.end();
        if (result.charge_number[i] != 0 && packing == bool(result.molecular_ion[i])) reject((packing ? "unconsumed component/" : "missing component/") + document.component_ids[i] + "/packing_diameter");
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
        for (const char* field : {"dipole_moment", "dipole_count", "quadrupole_moment", "quadrupole_count"}) {
            if ((field[0] == 'd' ? polar : quadrupolar) && component_fields.find({field, i}) == component_fields.end()) reject("missing component/" + document.component_ids[i] + "/" + field);
        }
    }
    for (std::size_t i = 0; i < components; ++i) for (std::size_t j = i + 1; j < components; ++j) {
        const bool same_sign = result.charge_number[i] * result.charge_number[j] > 0;
        const bool k_law = law_slots.count({1, i * components + j}) != 0;
        if (k_law && result.k_ij_reference_temperature_k[i * components + j] != 0.0) reject("unconsumed k_ij slope of a temperature law pair");
        if (same_sign && result.ionic_dispersion.is(IonicDispersion::ExcludeSameSign) && !(result.molecular_ion[i] && result.molecular_ion[j])) {
            if (k_law) reject("unconsumed k_ij temperature law under exclude-same-sign-ion-pairs");
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
        if (!k_law && result.k_ij[i * components + j] == 0.0 && result.k_ij[j * components + i] == 0.0) {
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
        const auto reverse = rule.index[1] * sites + rule.index[0];
        result.association_active[output] = result.association_active[reverse] = rule.identity == "wolbach-sandler" ? 2 : 1;
        // A declared correction or self-edge law evaluates this edge at each temperature (association_volume_correction).
        const auto first = canonical_slot(rule.index[2], rule.index[3], sites), second = canonical_slot(rule.index[4], rule.index[5], sites);
        const auto correction = corrections.find(canonical_slot(rule.index[0], rule.index[1], sites));
        if (correction != corrections.end() || law_slots.count({3, first}) || law_slots.count({4, first})
            || law_slots.count({3, second}) || law_slots.count({4, second})) {
            result.association_thermal_rule.push_back({static_cast<int>(output), static_cast<int>(reverse), static_cast<int>(first), static_cast<int>(second),
                correction == corrections.end() ? -1 : static_cast<int>(correction->second)});
            if (correction != corrections.end()) corrections.erase(correction);
            continue;
        }
        const double energy = 0.5 * (result.association_energy_over_k[self1] + result.association_energy_over_k[self2]);
        const double volume = std::sqrt(result.association_volume[self1] * result.association_volume[self2]);
        result.association_energy_over_k[output] = energy;
        result.association_energy_over_k[reverse] = energy;
        result.association_volume[output] = volume;
        result.association_volume[reverse] = volume;
    }
    if (!corrections.empty()) reject("unconsumed association_volume_correction");
    // Full EOS §5: omitted cross edges stay inactive, so each site-bearing component pair declares a cross
    // edge or no_cross_association; nothing is inferred.
    std::map<std::pair<std::size_t, std::size_t>, bool> cross;  // site-bearing component pair -> any active cross edge
    for (std::size_t s = 0; s < sites; ++s) {
        for (std::size_t t = 0; t < sites; ++t) {
            const std::pair<std::size_t, std::size_t> pair(document.site_component[s], document.site_component[t]);
            if (pair.first < pair.second) cross[pair] = cross[pair] || result.association_active[s * sites + t] != 0;
        }
    }
    for (const auto& [pair, active] : cross)
        if (no_cross.count(pair) == static_cast<std::size_t>(active))
            reject((active ? "no_cross_association contradicts a cross edge: " : "missing cross association or no_cross_association: ")
                   + document.component_ids[pair.first] + "/" + document.component_ids[pair.second]);
    // Selection fixes topology (full EOS §12): an uncharged selection has no
    // ionic term, so its packet-level electrolyte choices have no consumer.
    if (std::none_of(result.charge_number.begin(), result.charge_number.end(), [](int z) { return z != 0; })) {
        result.electrolyte_form = {};
        result.ionic_dispersion = {};
        result.permittivity_form = {};
        result.permittivity_solvent.clear();
        result.permittivity_solvent_normalized = false;
    }
    validate_cation_mixing(result);
    return result;
}

// NASA9 coefficients in units of the source gas constant, rescaled once to the Engine's.
ThermochemistryRecord prepare_thermochemistry(ThermochemistryRecord record, double source_gas_constant) {
    core::Scope scope("prepare_thermochemistry");
    if (!std::isfinite(source_gas_constant) || source_gas_constant <= 0.0) reject("invalid source gas constant");
    const double scale = source_gas_constant / GAS_CONSTANT_J_PER_MOL_K;
    for (auto& intervals : record.component_intervals) {
        validate_ideal_intervals(intervals);  // also refuses a species without intervals
        for (auto& interval : intervals) {
            auto* nasa9 = std::get_if<IdealNASA9>(&interval.correlation.form);
            if (nasa9 == nullptr) reject("thermochemistry packet interval is not NASA9");
            for (double& coefficient : nasa9->coefficients) coefficient *= scale;
        }
    }
    return record;
}

ParameterDocument fitted_document(const ParameterDocument& start, const Parameters& fitted) {
    core::Scope scope("fitted_document");
    static const std::array<std::string, detail::PARAMETER_FAMILY_COUNT> fields({"segment_count", "sigma_constant",
        "dispersion_energy_over_k", "k_ij", "l_ij", "association_energy_over_k", "association_volume", "born_diameter",
        "debye_huckel_diameter", "packing_diameter", "ion_suppression_coefficient", "born_c_dielectric",
        "ionic_region_relative_permittivity", "k_ij_reciprocal_temperature_slope", "dielectric_volume",
        "dielectric_temperature", "cation_pair_mixing", "cation_triple_mixing", "temperature_law_coefficient"});  // ParameterFamily order
    const Parameters base = build_parameters(start);
    const auto charged = std::find_if(base.charge_number.begin(), base.charge_number.end(), [](int z) { return z != 0; });
    ParameterDocument result = start;
    std::map<std::pair<std::size_t, std::size_t>, std::string> units;
    const auto& laws = base.temperature_law_coefficient_identity;  // a law term's slot is its identity's position
    for (auto& entry : result.entries) {
        std::size_t k = entry.index.empty() ? 0 : entry.index[0];
        if (entry.field == "cation_pair_mixing" || entry.field == "cation_triple_mixing") k = 0;
        if (entry.index.size() == 2)  // a site pair or a component pair
            k = k * (entry.field.rfind("association_", 0) == 0 ? start.site_component.size() : start.component_ids.size()) + entry.index[1];
        // The builder spreads the scalar suppression coefficient over the charged slots.
        if (entry.field == "ion_suppression_coefficient" && entry.index.empty()) {
            if (charged == base.charge_number.end()) continue;
            k = static_cast<std::size_t>(charged - base.charge_number.begin());
        }
        auto f = static_cast<std::size_t>(std::find(fields.begin(), fields.end(), entry.field) - fields.begin());
        const bool term = entry.field.find("temperature_law_") == 0 || entry.field == "association_volume_correction";
        const auto law = term ? std::find(laws.begin(), laws.end(), entry.identity) : laws.end();
        if (law != laws.end()) f = fields.size() - 1, k = static_cast<std::size_t>(law - laws.begin());
        if (f == fields.size()) continue;
        const auto family = static_cast<detail::ParameterFamily>(f);
        units[{f, k}] = entry.unit;
        const auto now = detail::family_values(fitted, family), was = detail::family_values(base, family);
        if (now.second != was.second || std::memcmp(now.first + k, was.first + k, sizeof(double)) == 0) continue;
        entry.magnitude = now.first[k] / (entry.unit == "liter / mole" ? unit_factor(entry, "molar volume") : 1.0);
    }
    const Parameters rebuilt = build_parameters(result);
    for (std::size_t f = 0; f < fields.size(); ++f) {
        const auto family = static_cast<detail::ParameterFamily>(f);
        const auto want = detail::family_values(fitted, family), got = detail::family_values(rebuilt, family);
        if (want.second != got.second) reject("fitted " + fields[f] + " slot count differs from the document's");
        for (std::size_t k = 0; k < want.second; ++k) {
            if (std::memcmp(want.first + k, got.first + k, sizeof(double)) == 0) continue;
            const auto unit = units.find({f, k});
            reject("fitted " + fields[f] + " slot " + std::to_string(k)
                + (unit == units.end() ? " has no document entry" : " does not round-trip in " + unit->second));
        }
    }
    return result;
}

}  // namespace epcsaft
