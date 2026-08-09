#include "eos.hpp"

#include "association.hpp"
#include "born.hpp"
#include "density.hpp"
#include "derivatives.hpp"
#include "electrolyte.hpp"
#include "electrolyte_reference.hpp"
#include "helmholtz.hpp"
#include "phase_block.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace epcsaft {
namespace {

// EqExclude: input and domain validation
bool positive_finite(const std::vector<double> &values) {
  return std::all_of(values.begin(), values.end(), [](double value) {
    return std::isfinite(value) && value > 0.0;
  });
}

// EqExclude: input and domain validation
// Explicit zero dispersion is a valid source parameter for ions.
bool nonnegative_finite(const std::vector<double> &values) {
  return std::all_of(values.begin(), values.end(), [](double value) {
    return std::isfinite(value) && value >= 0.0;
  });
}

// EqExclude: input and domain validation
bool finite(const std::vector<double> &values) {
  return std::all_of(values.begin(), values.end(),
                     [](double value) { return std::isfinite(value); });
}

std::vector<std::string> split_identity(const std::string &identity) {
  std::vector<std::string> parts;
  std::size_t begin = 0;
  while (begin <= identity.size()) {
    const std::size_t end = identity.find('/', begin);
    parts.push_back(identity.substr(begin, end - begin));
    if (end == std::string::npos) {
      break;
    }
    begin = end + 1;
  }
  return parts;
}

std::size_t named_index(const std::vector<std::string> &names,
                        const std::string &name, const char *kind) {
  const auto found = std::find(names.begin(), names.end(), name);
  if (found == names.end()) {
    throw std::invalid_argument(std::string("unknown ") + kind +
                                " in active parameter identity: " + name);
  }
  return static_cast<std::size_t>(found - names.begin());
}

// EqExclude: input and domain validation
bool polar_composition_admitted(const NativeParameters &parameters,
                                const std::vector<double> &mole_fractions) {
  if (parameters.polar_formulation == POLAR_FORMULATION_NONE) {
    return true;
  }
  for (std::size_t component = 0; component < mole_fractions.size();
       ++component) {
    if (mole_fractions[component] < 0.0) {
      return false;
    }
    const bool carries_moment =
        parameters.dipole_moment_debye[component] != 0.0 ||
        parameters.quadrupole_moment_debye_angstrom[component] != 0.0;
    if (carries_moment && mole_fractions[component] <= 0.0) {
      return false;
    }
  }
  return true;
}

std::optional<double> single_solvent_log_fugacity(
    const NativeModel &model, double temperature, double pressure,
    const std::vector<double> &solvent_fractions, const std::string &phase) {
  std::optional<std::size_t> solvent_index;
  for (std::size_t component = 0; component < solvent_fractions.size();
       ++component) {
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
  const auto state =
      model.evaluate_pressure(temperature, pressure, solvent_fractions, phase);
  if (!state.log_fugacity_coefficient.has_value()) {
    throw std::runtime_error("pure-solvent reference fugacity is unavailable");
  }
  return (*state.log_fugacity_coefficient)[*solvent_index];
}

} // namespace

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
      parameters_.packing_diameter_angstrom.size() != components ||
      parameters_.debye_huckel_diameter_angstrom.size() != components ||
      parameters_.born_diameter_angstrom.size() != components ||
      parameters_.solvation_factor.size() != components ||
      parameters_.dielectric_volume_m3_per_mol.size() != components ||
      parameters_.dielectric_temperature_k.size() != components ||
      parameters_.ion_specific_suppression_coefficient.size() != components ||
      parameters_.dipole_scaling_debye2.size() != components ||
      parameters_.polarizability_scaling_angstrom3.size() != components ||
      parameters_.correlation_integral_parameter.size() != components ||
      parameters_.dipole_moment_debye.size() != components ||
      parameters_.quadrupole_moment_debye_angstrom.size() != components) {
    throw std::invalid_argument("inconsistent native model dimensions");
  }
  if (!positive_finite(parameters_.segment_count) ||
      !positive_finite(parameters_.sigma_constant) ||
      !nonnegative_finite(parameters_.dispersion_energy_over_k) ||
      !finite(parameters_.k_ij) || !finite(parameters_.l_ij) ||
      !finite(parameters_.sigma_term_amplitude) ||
      !finite(parameters_.sigma_term_exponent) ||
      !finite(parameters_.association_energy_over_k) ||
      !finite(parameters_.association_volume) ||
      !finite(parameters_.molar_mass_kg_per_mol) ||
      !finite(parameters_.relative_permittivity) ||
      !std::isfinite(parameters_.ion_fraction_suppression_coefficient) ||
      !finite(parameters_.packing_diameter_angstrom) ||
      !finite(parameters_.debye_huckel_diameter_angstrom) ||
      !finite(parameters_.born_diameter_angstrom) ||
      !finite(parameters_.solvation_factor) ||
      !finite(parameters_.dielectric_volume_m3_per_mol) ||
      !finite(parameters_.dielectric_temperature_k) ||
      !finite(parameters_.ion_specific_suppression_coefficient) ||
      !finite(parameters_.dipole_scaling_debye2) ||
      !finite(parameters_.polarizability_scaling_angstrom3) ||
      !finite(parameters_.correlation_integral_parameter) ||
      !finite(parameters_.dipole_moment_debye) ||
      !finite(parameters_.quadrupole_moment_debye_angstrom) ||
      !std::isfinite(parameters_.ionic_region_relative_permittivity) ||
      !std::isfinite(parameters_.born_c_shell) ||
      !std::isfinite(parameters_.born_c_dielectric)) {
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
  const bool ionic = std::any_of(parameters_.charge_number.begin(),
                                 parameters_.charge_number.end(),
                                 [](int charge) { return charge != 0; });
  if (parameters_.polar_formulation != POLAR_FORMULATION_NONE &&
      parameters_.polar_formulation !=
          POLAR_FORMULATION_GROSS_VRABEC_POINT_MULTIPOLE) {
    throw std::invalid_argument("invalid polar formulation");
  }
  const bool has_multipole =
      std::any_of(parameters_.dipole_moment_debye.begin(),
                  parameters_.dipole_moment_debye.end(),
                  [](double value) { return value != 0.0; }) ||
      std::any_of(parameters_.quadrupole_moment_debye_angstrom.begin(),
                  parameters_.quadrupole_moment_debye_angstrom.end(),
                  [](double value) { return value != 0.0; });
  if (parameters_.polar_formulation == POLAR_FORMULATION_NONE &&
      has_multipole) {
    throw std::invalid_argument(
        "physical multipoles require a polar formulation");
  }
  for (std::size_t component = 0; component < components; ++component) {
    if (parameters_.dipole_moment_debye[component] < 0.0 ||
        parameters_.quadrupole_moment_debye_angstrom[component] < 0.0 ||
        (parameters_.charge_number[component] != 0 &&
         (parameters_.dipole_moment_debye[component] != 0.0 ||
          parameters_.quadrupole_moment_debye_angstrom[component] != 0.0))) {
      throw std::invalid_argument("invalid physical multipole parameters");
    }
  }
  if (ionic) {
    if (parameters_.electrolyte_model != ELECTROLYTE_MODEL_BORN &&
        parameters_.electrolyte_model !=
            ELECTROLYTE_MODEL_FULLY_DISSOCIATED_DEBYE_HUCKEL) {
      throw std::invalid_argument("invalid electrolyte model selection");
    }
    const bool born = parameters_.electrolyte_model == ELECTROLYTE_MODEL_BORN;
    const bool fully_dissociated =
        parameters_.electrolyte_model ==
        ELECTROLYTE_MODEL_FULLY_DISSOCIATED_DEBYE_HUCKEL;
    const bool valid_born_permittivity =
        parameters_.relative_permittivity_formulation ==
            RELATIVE_PERMITTIVITY_FORMULATION_ION_FRACTION_SUPPRESSION ||
        parameters_.relative_permittivity_formulation ==
            RELATIVE_PERMITTIVITY_FORMULATION_COMPONENT_MIXING ||
        parameters_.relative_permittivity_formulation ==
            RELATIVE_PERMITTIVITY_FORMULATION_DENSITY_VOLUME_TEMPERATURE ||
        parameters_.relative_permittivity_formulation ==
            RELATIVE_PERMITTIVITY_FORMULATION_ION_SPECIFIC_SUPPRESSION ||
        parameters_.relative_permittivity_formulation ==
            RELATIVE_PERMITTIVITY_FORMULATION_DIPOLAR_DENSITY_PERTURBATION;
    if ((born && !valid_born_permittivity) ||
        (fully_dissociated &&
         parameters_.relative_permittivity_formulation !=
             RELATIVE_PERMITTIVITY_FORMULATION_SOLVENT_ONLY)) {
      throw std::invalid_argument(
          "incompatible electrolyte and permittivity models");
    }
    std::size_t solvent_count = 0;
    for (std::size_t component = 0; component < components; ++component) {
      const bool solvent = parameters_.charge_number[component] == 0;
      solvent_count += static_cast<std::size_t>(solvent);
      const bool needs_tabulated_permittivity =
          parameters_.relative_permittivity_formulation ==
              RELATIVE_PERMITTIVITY_FORMULATION_ION_FRACTION_SUPPRESSION ||
          parameters_.relative_permittivity_formulation ==
              RELATIVE_PERMITTIVITY_FORMULATION_COMPONENT_MIXING ||
          parameters_.relative_permittivity_formulation ==
              RELATIVE_PERMITTIVITY_FORMULATION_ION_SPECIFIC_SUPPRESSION ||
          parameters_.relative_permittivity_formulation ==
              RELATIVE_PERMITTIVITY_FORMULATION_SOLVENT_ONLY;
      if (solvent &&
          (!(parameters_.molar_mass_kg_per_mol[component] > 0.0) ||
           (needs_tabulated_permittivity &&
            !(parameters_.relative_permittivity[component] > 0.0)))) {
        throw std::invalid_argument("invalid solvent dielectric parameters");
      }
      const double packing_diameter =
          parameters_.packing_diameter_angstrom[component];
      const double debye_diameter =
          parameters_.debye_huckel_diameter_angstrom[component];
      if ((solvent &&
           (packing_diameter != 0.0 || debye_diameter != 0.0)) ||
          (!solvent && !(packing_diameter > 0.0)) ||
          (!solvent && !(debye_diameter > 0.0))) {
        throw std::invalid_argument("invalid ionic diameter declaration");
      }
      if (born) {
        if (!(parameters_.ionic_region_relative_permittivity > 0.0) ||
            !(parameters_.solvation_factor[component] > 0.0)) {
          throw std::invalid_argument("invalid Born electrolyte parameters");
        }
        const double born_diameter = parameters_.born_diameter_angstrom[component];
        if ((solvent && born_diameter != 0.0) ||
            (!solvent && !(born_diameter > 0.0))) {
          throw std::invalid_argument("invalid Born diameter");
        }
      } else if (parameters_.electrolyte_model ==
                     ELECTROLYTE_MODEL_FULLY_DISSOCIATED_DEBYE_HUCKEL &&
                 ((parameters_.charge_number[component] != 0 &&
                  parameters_.relative_permittivity[component] != 0.0) ||
                  parameters_.born_diameter_angstrom[component] != 0.0 ||
                  parameters_.solvation_factor[component] != 0.0)) {
        throw std::invalid_argument(
            "invalid fully dissociated electrolyte parameters");
      }
      if (parameters_.relative_permittivity_formulation ==
              RELATIVE_PERMITTIVITY_FORMULATION_DENSITY_VOLUME_TEMPERATURE &&
          solvent &&
          (!(parameters_.dielectric_volume_m3_per_mol[component] > 0.0) ||
           !(parameters_.dielectric_temperature_k[component] > 0.0))) {
        throw std::invalid_argument(
            "invalid density formulation solvent parameters");
      }
      if (parameters_.relative_permittivity_formulation ==
              RELATIVE_PERMITTIVITY_FORMULATION_ION_SPECIFIC_SUPPRESSION &&
          !solvent &&
          !(parameters_.ion_specific_suppression_coefficient[component] >
            0.0)) {
        throw std::invalid_argument(
            "invalid ion-specific suppression parameter");
      }
      if (parameters_.relative_permittivity_formulation ==
          RELATIVE_PERMITTIVITY_FORMULATION_DIPOLAR_DENSITY_PERTURBATION) {
        if (parameters_.dipole_scaling_debye2[component] < 0.0 ||
            parameters_.polarizability_scaling_angstrom3[component] < 0.0 ||
            !(parameters_.correlation_integral_parameter[component] > 0.0)) {
          throw std::invalid_argument(
              "invalid dipolar perturbation parameters");
        }
      }
    }
    const bool valid_born_corrections =
        (parameters_.born_c_shell == 0.0 &&
         parameters_.born_c_dielectric == 0.0) ||
        (parameters_.born_c_shell == 1.0 &&
         parameters_.born_c_dielectric == 0.0) ||
        (parameters_.born_c_shell == 1.0 &&
         parameters_.born_c_dielectric == 1.0);
    if (born && !valid_born_corrections) {
      throw std::invalid_argument("invalid Born correction coefficients");
    }
    if (parameters_.relative_permittivity_formulation ==
            RELATIVE_PERMITTIVITY_FORMULATION_ION_FRACTION_SUPPRESSION &&
        !(parameters_.ion_fraction_suppression_coefficient > 0.0)) {
      throw std::invalid_argument("invalid ion-fraction suppression parameter");
    }
    if (parameters_.relative_permittivity_formulation ==
            RELATIVE_PERMITTIVITY_FORMULATION_ION_SPECIFIC_SUPPRESSION &&
        solvent_count != 1) {
      throw std::invalid_argument("mixed-solvent volume data are required");
    }
    if (fully_dissociated &&
        (parameters_.ion_fraction_suppression_coefficient != 0.0 ||
         parameters_.ionic_region_relative_permittivity != 0.0 ||
         parameters_.born_c_shell != 0.0 ||
         parameters_.born_c_dielectric != 0.0 ||
         std::any_of(parameters_.l_ij.begin(), parameters_.l_ij.end(),
                     [](double value) { return value != 0.0; }))) {
      throw std::invalid_argument(
          "invalid fully dissociated electrolyte parameters");
    }
  } else if (parameters_.electrolyte_model != ELECTROLYTE_MODEL_NONE ||
             parameters_.relative_permittivity_formulation !=
                 RELATIVE_PERMITTIVITY_FORMULATION_NONE) {
    throw std::invalid_argument("neutral model has electrolyte selection");
  }
}

NativeModel NativeModel::with_parameter_identities(
    const std::vector<std::string> &component_ids,
    const std::vector<std::string> &site_ids) const {
  if (component_ids.size() != parameters_.segment_count.size() ||
      site_ids.size() != parameters_.site_component.size() ||
      std::any_of(component_ids.begin(), component_ids.end(),
                  [](const std::string &value) { return value.empty(); }) ||
      std::any_of(site_ids.begin(), site_ids.end(),
                  [](const std::string &value) { return value.empty(); })) {
    throw std::invalid_argument("invalid native identity metadata");
  }
  for (std::size_t site = 0; site < site_ids.size(); ++site) {
    const auto parts = split_identity(site_ids[site]);
    const int component = parameters_.site_component[site];
    if (parts.size() != 2 || component < 0 ||
        static_cast<std::size_t>(component) >= component_ids.size() ||
        parts[0] != component_ids[static_cast<std::size_t>(component)] ||
        parts[1].empty()) {
      throw std::invalid_argument(
          "site identity metadata must use component/site identities");
    }
  }
  auto sorted_components = component_ids;
  auto sorted_sites = site_ids;
  std::sort(sorted_components.begin(), sorted_components.end());
  std::sort(sorted_sites.begin(), sorted_sites.end());
  if (std::adjacent_find(sorted_components.begin(), sorted_components.end()) !=
          sorted_components.end() ||
      std::adjacent_find(sorted_sites.begin(), sorted_sites.end()) !=
          sorted_sites.end()) {
    throw std::invalid_argument("duplicate native identity metadata");
  }
  NativeModel result(parameters_);
  result.component_ids_ = component_ids;
  result.site_ids_ = site_ids;
  return result;
}

// EqID: pure_active_parameter_mapping
NativeModel NativeModel::with_active_coordinates(
    const std::vector<NativeActiveCoordinate> &active_coordinates) const {
  if (active_coordinates.empty()) {
    throw std::invalid_argument("active parameter set must be nonempty");
  }
  NativeParameters updated = parameters_;
  for (const auto &coordinate : active_coordinates) {
    const auto parameter = resolve_active_coordinate(coordinate);
    switch (parameter.kind) {
    case NativeFixedTopologyParameterKind::SegmentCount:
      updated.segment_count[parameter.component] = parameter.value;
      break;
    case NativeFixedTopologyParameterKind::SegmentDiameter:
    case NativeFixedTopologyParameterKind::SegmentDiameterCorrelationConstant:
      updated.sigma_constant[parameter.component] = parameter.value;
      break;
    case NativeFixedTopologyParameterKind::SegmentDiameterCorrelationAmplitude:
      updated.sigma_term_amplitude[parameter.correlation_term] = parameter.value;
      break;
    case NativeFixedTopologyParameterKind::SegmentDiameterCorrelationExponent:
      updated.sigma_term_exponent[parameter.correlation_term] = parameter.value;
      break;
    case NativeFixedTopologyParameterKind::DispersionEnergyOverK:
      updated.dispersion_energy_over_k[parameter.component] = parameter.value;
      break;
    case NativeFixedTopologyParameterKind::AssociationEnergyOverK:
    case NativeFixedTopologyParameterKind::AssociationVolume: {
      const std::size_t sites = updated.site_component.size();
      const std::size_t forward = parameter.site_a * sites + parameter.site_b;
      const std::size_t reverse = parameter.site_b * sites + parameter.site_a;
      auto &values =
          parameter.kind == NativeFixedTopologyParameterKind::AssociationEnergyOverK
              ? updated.association_energy_over_k
              : updated.association_volume;
      values[forward] = parameter.value;
      values[reverse] = parameter.value;
      break;
    }
    case NativeFixedTopologyParameterKind::PairDispersion:
    case NativeFixedTopologyParameterKind::PairSize: {
      const std::size_t components = updated.segment_count.size();
      const std::size_t forward =
          parameter.component * components + parameter.component_b;
      const std::size_t reverse =
          parameter.component_b * components + parameter.component;
      auto &values =
          parameter.kind == NativeFixedTopologyParameterKind::PairDispersion
              ? updated.k_ij
              : updated.l_ij;
      values[forward] = parameter.value;
      values[reverse] = parameter.value;
      break;
    }
    case NativeFixedTopologyParameterKind::CrossAssociationEnergyCorrection: {
      const std::size_t sites = updated.site_component.size();
      for (std::size_t index = 0; index < parameter.association_edges.size(); ++index) {
        const auto edge = parameter.association_edges[index];
        const double energy = parameter.association_energy_scales[index] *
                              (1.0 - parameter.value);
        updated.association_energy_over_k[edge[0] * sites + edge[1]] = energy;
        updated.association_energy_over_k[edge[1] * sites + edge[0]] = energy;
      }
      break;
    }
    case NativeFixedTopologyParameterKind::RelativePermittivity:
      updated.relative_permittivity[parameter.component] = parameter.value;
      break;
    case NativeFixedTopologyParameterKind::BornDiameter:
      updated.born_diameter_angstrom[parameter.component] = parameter.value;
      break;
    case NativeFixedTopologyParameterKind::SolvationFactor:
    case NativeFixedTopologyParameterKind::SolvationFactorCorrelationConstant:
      updated.solvation_factor[parameter.component] = parameter.value;
      break;
    case NativeFixedTopologyParameterKind::IonFractionSuppressionCoefficient:
      updated.ion_fraction_suppression_coefficient = parameter.value;
      break;
    case NativeFixedTopologyParameterKind::IonicRegionRelativePermittivity:
      updated.ionic_region_relative_permittivity = parameter.value;
      break;
    }
  }
  NativeModel result(std::move(updated));
  result.component_ids_ = component_ids_;
  result.site_ids_ = site_ids_;
  return result;
}

NativeFixedTopologyParameter NativeModel::resolve_fixed_topology_parameter(
    const std::string &identity, double value) const {
  if (!std::isfinite(value)) {
    throw std::invalid_argument("active parameter values must be finite");
  }
  if (component_ids_.size() != parameters_.segment_count.size() ||
      site_ids_.size() != parameters_.site_component.size()) {
    throw std::logic_error(
        "canonical active parameters require native identity metadata");
  }
  const auto parts = split_identity(identity);
  if (parts.size() == 3 && parts[0] == "component") {
    const std::size_t component =
        named_index(component_ids_, parts[1], "component");
    NativeFixedTopologyParameterKind kind;
    if (parts[2] == "segment_count") {
      kind = NativeFixedTopologyParameterKind::SegmentCount;
    } else if (parts[2] == "segment_diameter") {
      kind = NativeFixedTopologyParameterKind::SegmentDiameter;
    } else if (parts[2] == "dispersion_energy_over_k") {
      kind = NativeFixedTopologyParameterKind::DispersionEnergyOverK;
    } else if (parts[2] == "relative_permittivity") {
      kind = NativeFixedTopologyParameterKind::RelativePermittivity;
    } else if (parts[2] == "born_diameter") {
      kind = NativeFixedTopologyParameterKind::BornDiameter;
    } else if (parts[2] == "solvation_factor") {
      kind = NativeFixedTopologyParameterKind::SolvationFactor;
    } else {
      throw std::logic_error("unsupported active parameter identity: " +
                             identity);
    }
    if (value <= 0.0) {
      throw std::invalid_argument(
          "direct component active values must be positive");
    }
    if ((kind == NativeFixedTopologyParameterKind::BornDiameter &&
         (parameters_.electrolyte_model != ELECTROLYTE_MODEL_BORN ||
          parameters_.charge_number[component] == 0)) ||
        (kind == NativeFixedTopologyParameterKind::SolvationFactor &&
         parameters_.electrolyte_model != ELECTROLYTE_MODEL_BORN) ||
        (kind == NativeFixedTopologyParameterKind::RelativePermittivity &&
         !(parameters_.relative_permittivity[component] > 0.0))) {
      throw std::invalid_argument(
          "component active coordinate is outside the model");
    }
    return {kind, component, 0, 0, value};
  }
  if ((parts.size() == 5 || parts.size() == 6) && parts[0] == "component" &&
      parts[2] == "segment_diameter") {
    const std::size_t component =
        named_index(component_ids_, parts[1], "component");
    if (parts.size() == 5 && parts[4] == "constant") {
      if (value <= 0.0) {
        throw std::invalid_argument(
            "segment-diameter correlation constant must be positive");
      }
      return {
          NativeFixedTopologyParameterKind::SegmentDiameterCorrelationConstant,
          component,
          0,
          0,
          value,
      };
    }
    if (parts.size() == 6 && parts[4].compare(0, 5, "term-") == 0 &&
        (parts[5] == "amplitude" || parts[5] == "exponent_coefficient")) {
      const std::string ordinal_text = parts[4].substr(5);
      if (ordinal_text.empty() ||
          !std::all_of(
              ordinal_text.begin(), ordinal_text.end(),
              [](char value) { return value >= '0' && value <= '9'; })) {
        throw std::logic_error(
            "invalid segment-diameter correlation term identity");
      }
      std::size_t ordinal = 0;
      for (const char value : ordinal_text) {
        const std::size_t digit = static_cast<std::size_t>(value - '0');
        if (ordinal > (std::numeric_limits<std::size_t>::max() - digit) / 10) {
          throw std::logic_error(
              "segment-diameter correlation term index is too large");
        }
        ordinal = ordinal * 10 + digit;
      }
      std::size_t observed = 0;
      std::size_t term = parameters_.sigma_term_component.size();
      for (std::size_t index = 0;
           index < parameters_.sigma_term_component.size(); ++index) {
        if (parameters_.sigma_term_component[index] !=
            static_cast<int>(component)) {
          continue;
        }
        if (observed == ordinal) {
          term = index;
          break;
        }
        ++observed;
      }
      if (term == parameters_.sigma_term_component.size()) {
        throw std::logic_error(
            "segment-diameter correlation term is outside the model");
      }
      const auto kind = parts[5] == "amplitude"
                            ? NativeFixedTopologyParameterKind::
                                  SegmentDiameterCorrelationAmplitude
                            : NativeFixedTopologyParameterKind::
                                  SegmentDiameterCorrelationExponent;
      return {kind, component, 0, 0, value, term};
    }
  }
  if (parts.size() == 6 && parts[0] == "association") {
    const std::size_t site_a = named_index(site_ids_, parts[1] + "/" + parts[2],
                                           "association endpoint");
    const std::size_t site_b = named_index(site_ids_, parts[3] + "/" + parts[4],
                                           "association endpoint");
    NativeFixedTopologyParameterKind kind;
    if (parts[5] == "energy_over_k") {
      kind = NativeFixedTopologyParameterKind::AssociationEnergyOverK;
    } else if (parts[5] == "volume") {
      kind = NativeFixedTopologyParameterKind::AssociationVolume;
    } else {
      throw std::logic_error("unsupported active parameter identity: " +
                             identity);
    }
    if (value <= 0.0) {
      throw std::invalid_argument("association active values must be positive");
    }
    return {
        kind,
        0,
        std::min(site_a, site_b),
        std::max(site_a, site_b),
        value,
    };
  }
  if (parts.size() == 5 && parts[0] == "component" &&
      parts[2] == "relative_permittivity" && parts[4] == "constant") {
    const std::size_t component =
        named_index(component_ids_, parts[1], "component");
    if (value <= 0.0 ||
        component >= parameters_.relative_permittivity.size() ||
        !(parameters_.relative_permittivity[component] > 0.0)) {
      throw std::invalid_argument(
          "relative-permittivity active value must be positive and applicable");
    }
    return {NativeFixedTopologyParameterKind::RelativePermittivity,
            component, 0, 0, value};
  }
  if (parts.size() == 5 && parts[0] == "component" &&
      parts[2] == "solvation_factor" && parts[4] == "constant") {
    const std::size_t component =
        named_index(component_ids_, parts[1], "component");
    if (value <= 0.0 || component >= parameters_.solvation_factor.size() ||
        parameters_.electrolyte_model != ELECTROLYTE_MODEL_BORN) {
      throw std::invalid_argument(
          "solvation-factor active value must be positive and applicable");
    }
    return {NativeFixedTopologyParameterKind::
                SolvationFactorCorrelationConstant,
            component, 0, 0, value};
  }
  if (identity ==
      "model/relative_permittivity/ion_fraction_suppression_coefficient") {
    if (value <= 0.0 ||
        parameters_.relative_permittivity_formulation != 1) {
      throw std::invalid_argument(
          "ion-fraction-suppression coordinate is outside the model");
    }
    return {NativeFixedTopologyParameterKind::IonFractionSuppressionCoefficient,
            0, 0, 0, value};
  }
  if (identity ==
      "model/relative_permittivity/ionic_region_relative_permittivity") {
    if (value <= 0.0 || parameters_.electrolyte_model != ELECTROLYTE_MODEL_BORN) {
      throw std::invalid_argument(
          "ionic-region-permittivity coordinate is outside the model");
    }
    return {NativeFixedTopologyParameterKind::IonicRegionRelativePermittivity,
            0, 0, 0, value};
  }
  throw std::logic_error("unsupported active parameter identity: " +
                         identity);
}

NativeFixedTopologyParameter NativeModel::resolve_active_coordinate(
    const NativeActiveCoordinate &coordinate) const {
  if (!std::isfinite(coordinate.value)) {
    throw std::invalid_argument("active coordinate value must be finite");
  }
  if (coordinate.application == "reciprocal-cross-dispersion" ||
      coordinate.application == "reciprocal-cross-size" ||
      coordinate.application == "role-aware-cross-association-energy") {
    if (coordinate.indices.size() != 2 ||
        coordinate.indices[0] == coordinate.indices[1] ||
        coordinate.indices[0] >= component_ids_.size() ||
        coordinate.indices[1] >= component_ids_.size()) {
      throw std::invalid_argument(
          "pair coordinate endpoints are not canonical model indices");
    }
    const auto parts = split_identity(coordinate.identity);
    const auto endpoint_ids = std::minmax(component_ids_[coordinate.indices[0]],
                                          component_ids_[coordinate.indices[1]]);
    const auto endpoint_indices =
        std::minmax(coordinate.indices[0], coordinate.indices[1]);
    const std::string prefix =
        "pair/" + endpoint_ids.first + "/" + endpoint_ids.second + "/";
    NativeFixedTopologyParameterKind kind;
    if (coordinate.application == "reciprocal-cross-dispersion") {
      kind = NativeFixedTopologyParameterKind::PairDispersion;
      if (coordinate.family != "k_ij" || coordinate.identity != prefix + "k_ij") {
        throw std::invalid_argument("pair-dispersion descriptor is inconsistent");
      }
    } else if (coordinate.application == "reciprocal-cross-size") {
      kind = NativeFixedTopologyParameterKind::PairSize;
      if (coordinate.family != "l_ij" || coordinate.identity != prefix + "l_ij") {
        throw std::invalid_argument("pair-size descriptor is inconsistent");
      }
    } else {
      kind = NativeFixedTopologyParameterKind::CrossAssociationEnergyCorrection;
      if (coordinate.family != "k_hb_ij" || parts.size() != 7 ||
          parts[0] != "pair" || parts[1] != endpoint_ids.first ||
          parts[2] != endpoint_ids.second ||
          parts[3] != "association" || parts[4] > parts[5] ||
          parts[6] != "k_hb_ij" || coordinate.association_edges.empty() ||
          coordinate.association_edges.size() !=
              coordinate.association_energy_scales.size()) {
        throw std::invalid_argument(
            "cross-association descriptor is inconsistent");
      }
    }
    return {kind,
            endpoint_indices.first,
            0,
            0,
            coordinate.value,
            0,
            endpoint_indices.second,
            coordinate.association_edges,
            coordinate.association_energy_scales};
  }

  const auto resolved =
      resolve_fixed_topology_parameter(coordinate.identity, coordinate.value);
  std::string expected_application;
  std::vector<std::size_t> expected_indices;
  switch (resolved.kind) {
  case NativeFixedTopologyParameterKind::SegmentCount:
  case NativeFixedTopologyParameterKind::SegmentDiameter:
  case NativeFixedTopologyParameterKind::DispersionEnergyOverK:
  case NativeFixedTopologyParameterKind::BornDiameter:
  case NativeFixedTopologyParameterKind::SolvationFactor:
    expected_application = "component-property";
    expected_indices = {resolved.component};
    break;
  case NativeFixedTopologyParameterKind::SegmentDiameterCorrelationConstant:
  case NativeFixedTopologyParameterKind::SegmentDiameterCorrelationAmplitude:
  case NativeFixedTopologyParameterKind::SegmentDiameterCorrelationExponent:
  case NativeFixedTopologyParameterKind::SolvationFactorCorrelationConstant:
    expected_application = "typed-correlation-coefficient";
    expected_indices = {resolved.component};
    break;
  case NativeFixedTopologyParameterKind::RelativePermittivity:
    expected_application =
        split_identity(coordinate.identity).size() == 3
            ? "component-property"
            : "typed-correlation-coefficient";
    expected_indices = {resolved.component};
    break;
  case NativeFixedTopologyParameterKind::AssociationEnergyOverK:
  case NativeFixedTopologyParameterKind::AssociationVolume:
    expected_application = "symmetric-association-property";
    expected_indices = {resolved.site_a, resolved.site_b};
    break;
  case NativeFixedTopologyParameterKind::IonFractionSuppressionCoefficient:
  case NativeFixedTopologyParameterKind::IonicRegionRelativePermittivity:
    expected_application = "auxiliary-model-coefficient";
    break;
  case NativeFixedTopologyParameterKind::PairDispersion:
  case NativeFixedTopologyParameterKind::PairSize:
  case NativeFixedTopologyParameterKind::CrossAssociationEnergyCorrection:
    throw std::logic_error("pair descriptor did not use its resolved application");
  }
  if (coordinate.application != expected_application ||
      coordinate.indices != expected_indices ||
      !coordinate.association_edges.empty() ||
      !coordinate.association_energy_scales.empty()) {
    throw std::invalid_argument("active coordinate descriptor is inconsistent");
  }
  return resolved;
}

bool NativeModel::has_temperature_independent_segment_diameter(
    std::size_t component) const {
  if (component >= parameters_.sigma_constant.size()) {
    throw std::out_of_range("segment-diameter component is out of range");
  }
  return std::none_of(parameters_.sigma_term_component.begin(),
                      parameters_.sigma_term_component.end(),
                      [component](int term_component) {
                        return term_component == static_cast<int>(component);
                      });
}

// EqExclude: input and domain validation
bool NativeModel::supports_pure_transport() const {
  return detail::supports_pure_phase_block(parameters_);
}

// EqExclude: input and domain validation
bool NativeModel::supports_fixed_topology_association() const {
  return detail::supports_fixed_topology_associating_phase_block(parameters_);
}

// EqExclude: input and domain validation
bool NativeModel::supports_fixed_topology_associating_mixture()
    const {
  return detail::supports_fixed_topology_associating_mixture_phase_block(
      parameters_);
}

// EqExclude: input and domain validation
bool NativeModel::supports_mixture_transport() const {
  return detail::supports_mixture_phase_block(parameters_);
}

// EqExclude: input and domain validation
bool NativeModel::supports_electrolyte_transport() const {
  return detail::supports_electrolyte_phase_block(parameters_);
}

// EqExclude: serialization
std::size_t NativeModel::component_count() const {
  return parameters_.segment_count.size();
}

// EqExclude: serialization
const std::vector<int> &NativeModel::charge_numbers() const {
  return parameters_.charge_number;
}

// EqExclude: serialization
const std::vector<int> &NativeModel::site_components() const {
  return parameters_.site_component;
}

// EqExclude: serialization
const std::vector<int> &NativeModel::site_multiplicities() const {
  return parameters_.site_multiplicity;
}

// EqExclude: serialization
const std::vector<int> &NativeModel::association_active_pairs() const {
  return parameters_.association_active;
}

// EqExclude: serialization
const std::vector<double> &NativeModel::molar_masses_kg_per_mol() const {
  return parameters_.molar_mass_kg_per_mol;
}

// EqExclude: transport and buffer mechanics
NativePurePhaseBlock NativeModel::evaluate_pure_phase(double temperature,
                                                      double amount,
                                                      double volume) const {
  return detail::pure_phase_block(temperature, amount, volume, parameters_);
}

// EqExclude: transport and buffer mechanics
NativeFixedTopologyPhaseBlock NativeModel::evaluate_fixed_topology_association(
    double temperature, double amount, double volume,
    const std::vector<NativeFixedTopologyParameter> &active_parameters) const {
  return detail::fixed_topology_parameterized_pure_phase_block(
      temperature, amount, volume, active_parameters, parameters_);
}

// EqExclude: transport and buffer mechanics
NativeMixturePhaseBlock
NativeModel::evaluate_mixture_phase(double temperature,
                                    const std::vector<double> &amounts,
                                    double volume) const {
  if (supports_fixed_topology_associating_mixture()) {
    return detail::fixed_topology_parameterized_mixture_phase_block(
        temperature, amounts, volume, {}, parameters_);
  }
  return detail::mixture_phase_block(temperature, amounts, volume, parameters_);
}

// EqExclude: transport and buffer mechanics
NativeMixturePhaseBlock
NativeModel::evaluate_phase_derivatives(
    double temperature, const std::vector<double> &amounts, double volume,
    const std::vector<NativeActiveCoordinate> &coordinates,
    std::size_t derivative_order) const {
  if (derivative_order != 2) {
    throw std::invalid_argument(
        "phase derivative order must be exactly two");
  }
  if (component_ids_.size() != parameters_.segment_count.size() ||
      coordinates.empty()) {
    throw std::invalid_argument(
        "phase derivatives require resolved active coordinates");
  }
  std::vector<NativeFixedTopologyParameter> resolved;
  resolved.reserve(coordinates.size());
  for (const auto &coordinate : coordinates) {
    resolved.push_back(resolve_active_coordinate(coordinate));
  }
  return detail::fixed_topology_parameterized_mixture_phase_block(
      temperature, amounts, volume, resolved, parameters_);
}

// EqExclude: transport and buffer mechanics
NativeMixturePhaseBlock NativeModel::evaluate_electrolyte_phase(
    double temperature, const std::vector<double> &amounts, double volume,
    double ion_mole_fraction_max) const {
  return detail::electrolyte_phase_block(temperature, amounts, volume,
                                         ion_mole_fraction_max, parameters_);
}

// EqExclude: transport and buffer mechanics
NativeMixturePhaseBlock NativeModel::evaluate_parameterized_phase(
    double temperature, const std::vector<double> &amounts, double volume,
    double ion_mole_fraction_max,
    const std::vector<NativeActiveCoordinate> &active_parameters) const {
  if (supports_electrolyte_transport()) {
    detail::validate_electrolyte_phase_state(
        temperature, amounts, volume, ion_mole_fraction_max, parameters_);
  } else {
    detail::validate_mixture_phase_state(temperature, amounts, volume,
                                         parameters_);
  }
  return evaluate_phase_derivatives(temperature, amounts, volume,
                                    active_parameters, 2);
}

// EqExclude: transport and buffer mechanics
double NativeModel::evaluate_electrolyte_phase_value(
    double temperature, const std::vector<double> &amounts, double volume,
    double ion_mole_fraction_max) const {
  return evaluate_electrolyte_phase(temperature, amounts, volume,
                                    ion_mole_fraction_max)
      .helmholtz_over_rt_reference_amount;
}

// EqExclude: transport and buffer mechanics
NativePackingFractionBlock NativeModel::evaluate_packing_fraction(
    double temperature, const std::vector<double> &amounts, double volume,
    double ion_mole_fraction_max) const {
  return detail::packing_fraction_block(temperature, amounts, volume,
                                        ion_mole_fraction_max, parameters_);
}

// EqExclude: transport and buffer mechanics
NativePackingFractionBlock NativeModel::evaluate_packing_fraction_parameters(
    double temperature, const std::vector<double> &amounts, double volume,
    double ion_mole_fraction_max,
    const std::vector<NativeActiveCoordinate> &active_parameters) const {
  std::vector<NativeFixedTopologyParameter> resolved;
  resolved.reserve(active_parameters.size());
  for (const auto &coordinate : active_parameters) {
    resolved.push_back(resolve_active_coordinate(coordinate));
  }
  return detail::packing_fraction_block_active_parameters(
      temperature, amounts, volume, ion_mole_fraction_max, parameters_,
      resolved);
}

// EqExclude: transport and buffer mechanics
NativeInversePackingGeometryBlock
NativeModel::evaluate_inverse_packing_geometry(
    double temperature, const std::vector<double> &amounts,
    double log_packing_fraction,
    const std::vector<NativeActiveCoordinate> &active_parameters) const {
  std::vector<NativeFixedTopologyParameter> resolved;
  resolved.reserve(active_parameters.size());
  for (const auto &coordinate : active_parameters) {
    resolved.push_back(resolve_active_coordinate(coordinate));
  }
  return detail::inverse_packing_geometry_block(temperature, amounts,
                                                log_packing_fraction,
                                                parameters_, resolved);
}

// EqID: packing_volume_bounds
NativeMolarVolumeBoundsParameterBlock
NativeModel::molar_volume_bounds_parameters(
    double temperature, const std::vector<double> &mole_fractions,
    double packing_fraction_min, double packing_fraction_max,
    double ion_mole_fraction_max,
    const std::vector<NativeActiveCoordinate> &active_parameters) const {
  const auto bounds = with_active_coordinates(active_parameters)
                          .molar_volume_bounds(temperature, mole_fractions,
                                               packing_fraction_min,
                                               packing_fraction_max,
                                               ion_mole_fraction_max);
  const auto coefficient = evaluate_packing_fraction_parameters(
      temperature, mole_fractions, 1.0, ion_mole_fraction_max,
      active_parameters);
  const std::size_t parameter_offset = mole_fractions.size() + 1;
  std::vector<double> derivatives(2 * active_parameters.size());
  for (std::size_t parameter = 0; parameter < active_parameters.size();
       ++parameter) {
    const double coefficient_derivative =
        coefficient.gradient[parameter_offset + parameter];
    derivatives[parameter] = coefficient_derivative / packing_fraction_max;
    derivatives[active_parameters.size() + parameter] =
        coefficient_derivative / packing_fraction_min;
  }
  return {bounds, std::move(derivatives)};
}

// EqID: formula_unit_composition
std::vector<double> NativeModel::reference_composition(
    const std::vector<double> &solvent_fractions,
    const std::array<std::size_t, 2> &ion_indices,
    const std::array<int, 2> &stoichiometry, double molality) const {
  return detail::formula_composition(solvent_fractions,
                                     parameters_.molar_mass_kg_per_mol,
                                     ion_indices, stoichiometry, molality);
}

// EqID: reference_molality_sequence
NativeReferenceSequence NativeModel::electrolyte_reference_sequence(
    double temperature, double pressure,
    const std::vector<double> &solvent_fractions,
    const std::array<std::size_t, 2> &ion_indices,
    const std::array<int, 2> &stoichiometry, const std::string &phase,
    NativeEvaluationBudget *budget) const {
  auto checkpoint = [budget] {
    if (budget != nullptr) {
      budget->checkpoint();
    }
  };
  const auto result = detail::electrolyte_reference_sequence(
      solvent_fractions, parameters_.molar_mass_kg_per_mol, ion_indices,
      stoichiometry,
      [this, temperature, pressure, &ion_indices,
       &phase](const std::vector<double> &composition) {
        const auto state =
            evaluate_pressure(temperature, pressure, composition, phase);
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
      checkpoint);
  const auto pure_solvent_log_fugacity = single_solvent_log_fugacity(
      *this, temperature, pressure, solvent_fractions, phase);
  return {
      result.log_fugacity,
      pure_solvent_log_fugacity,
      result.molality,
      result.convergence_error,
  };
}

// EqID: packing_volume_bounds
std::array<double, 2> NativeModel::molar_volume_bounds(
    double temperature, const std::vector<double> &mole_fractions,
    double packing_fraction_min, double packing_fraction_max,
    double ion_mole_fraction_max) const {
  const std::size_t components = parameters_.segment_count.size();
  if (!std::isfinite(temperature) || temperature <= 0.0 ||
      mole_fractions.size() != components ||
      !std::isfinite(packing_fraction_min) ||
      !std::isfinite(packing_fraction_max) || packing_fraction_min <= 0.0 ||
      packing_fraction_min >= packing_fraction_max) {
    throw std::invalid_argument("invalid molar-volume bound inputs");
  }
  if (packing_fraction_max >= detail::close_packing_fraction()) {
    throw std::invalid_argument(
        "packing-fraction upper bound is outside the model domain");
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
  if (std::abs(total_fraction - 1.0) > 1.0e-12 ||
      std::abs(charge_fraction) > 1.0e-12 ||
      (std::isfinite(ion_mole_fraction_max) &&
       ion_fraction > ion_mole_fraction_max + 1.0e-12)) {
    throw std::invalid_argument("composition is outside the model domain");
  }

  const double coefficient = detail::packing_fraction_coefficient(
      temperature, mole_fractions, parameters_);
  if (!(coefficient > 0.0) || !std::isfinite(coefficient)) {
    throw std::invalid_argument("invalid model packing coefficient");
  }
  const std::array<double, 2> bounds = {
      coefficient / packing_fraction_max,
      coefficient / packing_fraction_min,
  };
  if (!(bounds[0] > 0.0) || !(bounds[0] < bounds[1]) ||
      !std::isfinite(bounds[1])) {
    throw std::invalid_argument("invalid model molar-volume bounds");
  }
  return bounds;
}

double NativeModel::molar_density_upper_bound(
    double temperature, const std::vector<double> &mole_fractions) const {
  const auto bounds =
      molar_volume_bounds(temperature, mole_fractions, 1.0e-12,
                          std::nextafter(detail::close_packing_fraction(), 0.0),
                          std::numeric_limits<double>::quiet_NaN());
  return 1.0 / bounds[0];
}

NativeResult
NativeModel::evaluate(double temperature, double molar_density,
                      const std::vector<double> &mole_fractions) const {
  const std::size_t components = parameters_.segment_count.size();
  if (!std::isfinite(temperature) || temperature <= 0.0 ||
      !std::isfinite(molar_density) || molar_density <= 0.0 ||
      mole_fractions.size() != components || !finite(mole_fractions)) {
    throw std::invalid_argument("invalid native state");
  }
  if (!polar_composition_admitted(parameters_, mole_fractions)) {
    throw std::invalid_argument("polar model mole fractions must be positive");
  }
  const auto sigma = detail::segment_sigma(temperature, parameters_);
  if (!positive_finite(sigma)) {
    throw std::invalid_argument("invalid segment diameter correlation");
  }
  const auto diameter =
      detail::segment_diameter(temperature, parameters_, sigma);
  const auto site_fraction = detail::solve_site_fractions(
      temperature, molar_density, mole_fractions, parameters_, sigma, diameter);
  const auto electrolyte = detail::electrolyte_contribution(
      temperature, molar_density, temperature, molar_density, mole_fractions,
      mole_fractions, parameters_);
  const double born = detail::born_contribution(
      temperature, mole_fractions, parameters_,
      electrolyte.bulk_relative_permittivity);
  const auto values = detail::residual_helmholtz(
      temperature, molar_density, mole_fractions, parameters_, sigma, diameter,
      site_fraction, electrolyte.debye_huckel, born);

  const double residual = values.hard_chain + values.dispersion +
                          values.association + values.polar +
                          values.debye_huckel + values.born;
  const auto derivatives = detail::thermodynamic_derivatives(
      temperature, molar_density, mole_fractions, parameters_, site_fraction);
  // EqID: compressibility_factor
  const double compressibility = 1.0 + molar_density * derivatives.density;
  double weighted_composition_derivative = 0.0;
  for (std::size_t component = 0; component < components; ++component) {
    weighted_composition_derivative +=
        mole_fractions[component] * derivatives.composition[component];
  }
  std::vector<double> chemical_potential(components);
  std::optional<std::vector<double>> log_fugacity;
  if (std::isfinite(compressibility) && compressibility > 0.0) {
    log_fugacity.emplace(components);
  }
  for (std::size_t component = 0; component < components; ++component) {
    // EqID: residual_chemical_potential
    chemical_potential[component] = residual + compressibility - 1.0 +
                                    derivatives.composition[component] -
                                    weighted_composition_derivative;
    if (log_fugacity.has_value()) {
      // EqID: log_fugacity_coefficient
      (*log_fugacity)[component] =
          chemical_potential[component] - std::log(compressibility);
    }
  }
  const double z_temperature = molar_density * derivatives.temperature_density;
  const double z_density =
      derivatives.density + molar_density * derivatives.density_density;
  const double pressure_slope_factor =
      compressibility + molar_density * z_density;
  const double cv_over_r =
      -2.0 * temperature * derivatives.temperature -
      temperature * temperature * derivatives.temperature_temperature;
  std::optional<double> residual_cp;
  if (std::isfinite(pressure_slope_factor) && pressure_slope_factor > 1.0e-10) {
    // EqID: residual_heat_capacity_cp
    const double cp_over_r =
        cv_over_r +
        (compressibility + temperature * z_temperature) *
            (compressibility + temperature * z_temperature) /
            pressure_slope_factor -
        1.0;
    if (std::isfinite(cp_over_r)) {
      residual_cp = detail::GAS_CONSTANT * cp_over_r;
    }
  }
  return {
      values.hard_chain,
      values.dispersion,
      values.association,
      values.polar,
      values.debye_huckel,
      values.born,
      electrolyte.bulk_relative_permittivity,
      residual,
      compressibility,
      // EqID: pressure
      compressibility * molar_density * detail::GAS_CONSTANT * temperature,
      std::move(chemical_potential),
      std::move(log_fugacity),
      // EqID: residual_entropy
      detail::GAS_CONSTANT *
          (-residual - temperature * derivatives.temperature),
      // EqID: residual_internal_energy
      -detail::GAS_CONSTANT * temperature * temperature *
          derivatives.temperature,
      // EqID: residual_enthalpy
      detail::GAS_CONSTANT * temperature *
          (-temperature * derivatives.temperature + compressibility - 1.0),
      // EqID: residual_heat_capacity_cv
      detail::GAS_CONSTANT * cv_over_r,
      residual_cp,
      molar_density,
      // EqID: pressure_density_derivative
      detail::GAS_CONSTANT * temperature * pressure_slope_factor,
      std::nullopt,
  };
}

NativeResult
NativeModel::evaluate_pressure(double temperature, double pressure,
                               const std::vector<double> &mole_fractions,
                               const std::string &phase) const {
  const std::size_t components = parameters_.segment_count.size();
  if (!std::isfinite(temperature) || temperature <= 0.0 ||
      !std::isfinite(pressure) || pressure <= 0.0 ||
      mole_fractions.size() != components || !finite(mole_fractions)) {
    throw std::invalid_argument("invalid native pressure state");
  }
  if (!polar_composition_admitted(parameters_, mole_fractions)) {
    throw std::invalid_argument("polar model mole fractions must be positive");
  }
  // EqID: density_upper_bound
  const double packing_coefficient = detail::packing_fraction_coefficient(
      temperature, mole_fractions, parameters_);
  if (!(packing_coefficient > 0.0) || !std::isfinite(packing_coefficient)) {
    throw std::invalid_argument("invalid pressure-state packing coefficient");
  }
  const auto root = detail::solve_density(
      pressure,
      // EqID: ideal_gas_density
      pressure / (detail::GAS_CONSTANT * temperature),
      detail::close_packing_fraction() / packing_coefficient, phase,
      [this, temperature, &mole_fractions](double density) {
        return evaluate(temperature, density, mole_fractions);
      });
  auto result = evaluate(temperature, root.molar_density, mole_fractions);
  result.density_diagnostics = NativeDensityDiagnostics{
      root.bracket_lower,
      root.bracket_upper,
      root.pressure_residual,
      root.pressure_density_derivative,
      root.iterations,
      root.stable,
      root.branch,
      root.stable_root_count,
      root.selected_stable_root_index,
  };
  return result;
}

} // namespace epcsaft
