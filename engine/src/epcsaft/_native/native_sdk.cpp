#include "native_sdk.hpp"
#include "density.hpp"
#include "phase_block.hpp"

#include <epcsaft/native_sdk_v1.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace epcsaft {
namespace {

// EqExclude: input and domain validation
bool active_born_model_admitted(
    const std::string& fingerprint,
    const std::vector<std::string>& component_ids,
    const std::vector<int>& charges,
    bool staged_trial_bundle_admitted
) {
    if (component_ids.size() != 3 || charges.size() != 3
        || component_ids[0] != "water" || charges[0] != 0) {
        return false;
    }
    struct Admission final {
        const char* fingerprint;
        const char* active_ion;
        const char* counter_ion;
        int active_charge;
    };
    constexpr Admission admissions[] = {
        {
            "sha256:dbed1c05a92d922369b9b631e89fe10add759e7c8d004ea53bc0a0c90d2d6f4a",
            "lithium-cation",
            "chloride-anion",
            1,
        },
        {
            "sha256:d01bc3e64f19293902c00d9cc52060fccc3945ce736f0cd3e2f756be226cf5bf",
            "sodium-cation",
            "chloride-anion",
            1,
        },
        {
            "sha256:05f2db21ad51a27e7b7b8fca267e8a38eb663145c980eae6bee7dce02700af5b",
            "potassium-cation",
            "chloride-anion",
            1,
        },
        {
            "sha256:a8383fb8f70a1bb8c258f1662c792cd5507af1d23db71ad4575026b5dcab55ea",
            "chloride-anion",
            "sodium-cation",
            -1,
        },
        {
            "sha256:f313c30e9471daed0638ab70ff43eb033ed828724547c4b32857442ed35f1107",
            "bromide-anion",
            "sodium-cation",
            -1,
        },
    };
    return std::any_of(
        std::begin(admissions),
        std::end(admissions),
        [&](const Admission& admission) {
            return component_ids[1] == admission.active_ion
                && component_ids[2] == admission.counter_ion
                && charges[1] == admission.active_charge
                && charges[2] == -admission.active_charge
                && (staged_trial_bundle_admitted
                    || fingerprint == admission.fingerprint);
        }
    );
}

// EqExclude: input and domain validation
bool aqueous_kij_model_admitted(
    const std::string& fingerprint,
    const std::vector<std::string>& component_ids,
    const std::vector<int>& charges,
    bool staged_trial_bundle_admitted
) {
    if (component_ids.size() != 3 || charges != std::vector<int>({0, 1, -1})
        || component_ids[0] != "water") {
        return false;
    }
    struct Admission final {
        const char* fingerprint;
        const char* cation;
        const char* anion;
    };
    constexpr Admission admissions[] = {
        {
            "sha256:dbed1c05a92d922369b9b631e89fe10add759e7c8d004ea53bc0a0c90d2d6f4a",
            "lithium-cation", "chloride-anion",
        },
        {
            "sha256:d01bc3e64f19293902c00d9cc52060fccc3945ce736f0cd3e2f756be226cf5bf",
            "sodium-cation", "chloride-anion",
        },
        {
            "sha256:05f2db21ad51a27e7b7b8fca267e8a38eb663145c980eae6bee7dce02700af5b",
            "potassium-cation", "chloride-anion",
        },
        {
            "sha256:19be2cea6620c8a6cd3c965325422bded799d2b39acb240a4c4e858059d6412e",
            "lithium-cation", "bromide-anion",
        },
        {
            "sha256:80f83701040d7d0863c7121b7379b86a9e8bd30c1fa3e4a5651330c41c588131",
            "sodium-cation", "bromide-anion",
        },
        {
            "sha256:8ec6c31ada39713132b340a184e0c06b0f9d070892be2aaa721e9aa92015ae31",
            "potassium-cation", "bromide-anion",
        },
    };
    return std::any_of(
        std::begin(admissions), std::end(admissions), [&](const Admission& item) {
            return component_ids[1] == item.cation
                && component_ids[2] == item.anion
                && (staged_trial_bundle_admitted
                    || fingerprint == item.fingerprint);
        }
    );
}

// EqExclude: input and domain validation
bool ion_fraction_suppression_model_admitted(
    const std::vector<std::string>& component_ids,
    const std::vector<int>& charges
) {
    return component_ids.size() == 3
        && charges == std::vector<int>({0, 1, -1})
        && (component_ids[0] == "water" || component_ids[0] == "methanol");
}

// EqExclude: input and domain validation
bool ion_solvation_kij_model_admitted(
    const std::string& parameter_fingerprint,
    const std::string& topology_fingerprint,
    const std::vector<std::string>& component_ids,
    const std::vector<int>& charges
) {
    struct Admission final {
        std::array<const char*, 3> component_ids;
        const char* parameter_fingerprint;
        const char* topology_fingerprint;
    };
    constexpr std::array<Admission, 2> admissions{{
        {
            {"methanol", "potassium-cation", "bromide-anion"},
            "sha256:57920a78527ba21924a8f3e3916e40a8b842fae5127c4cfc6ff3e1fade418699",
            "sha256:a6879554a81191d2210ece2eed240f251a0f45327fa3877978228362bcf80443",
        },
        {
            {"ethanol", "sodium-cation", "chloride-anion"},
            "sha256:02d7e56ffd746f91632c14c9e68c247df454b0db53514499dd26fbeedfb11016",
            "sha256:7998ddbc59d7059e7ac3ee0810d1b97d1d564264c98dea41b55c6c8192fe2a71",
        },
    }};
    return charges == std::vector<int>({0, 1, -1})
        && std::any_of(
            admissions.begin(),
            admissions.end(),
            [&](const Admission& admission) {
                return component_ids
                        == std::vector<std::string>(
                            admission.component_ids.begin(),
                            admission.component_ids.end()
                        )
                    && parameter_fingerprint
                        == admission.parameter_fingerprint
                    && topology_fingerprint
                        == admission.topology_fingerprint;
            }
        );
}

// EqExclude: input and domain validation
bool aqueous_solvation_factor_model_admitted(
    const std::string& fingerprint,
    const std::vector<std::string>& component_ids,
    const std::vector<int>& charges,
    bool staged_trial_bundle_admitted
) {
    constexpr const char* sodium_bromide_fingerprint =
        "sha256:80f83701040d7d0863c7121b7379b86a9e8bd30c1fa3e4a5651330c41c588131";
    return component_ids == std::vector<std::string>({
               "water", "sodium-cation", "bromide-anion"
           })
        && charges == std::vector<int>({0, 1, -1})
        && (staged_trial_bundle_admitted
            || fingerprint == sodium_bromide_fingerprint);
}

// EqExclude: input and domain validation
bool neutral_reference_model_shape_admitted(
    const std::vector<std::string>& component_ids,
    const std::vector<int>& charges
) {
    if (component_ids.size() != charges.size()) {
        return false;
    }
    std::size_t water_count = 0;
    std::size_t cation_count = 0;
    std::size_t anion_count = 0;
    for (std::size_t component = 0; component < component_ids.size(); ++component) {
        water_count += component_ids[component] == "water" && charges[component] == 0;
        cation_count += charges[component] > 0;
        anion_count += charges[component] < 0;
    }
    return water_count == 1 && cation_count > 0 && anion_count > 0;
}

struct NativeSdkEndpoint final {
    NativeModel model;
    std::string parameter_fingerprint;
    std::string topology_fingerprint;
    double temperature_min_k;
    double temperature_max_k;
    double pressure_min_pa;
    double pressure_max_pa;
    double ion_mole_fraction_max;
    bool phase_callbacks_admitted;
    bool active_born_admitted;
    bool aqueous_kij_admitted;
    bool ion_fraction_suppression_admitted;
    bool ion_solvation_kij_admitted;
    bool aqueous_solvation_factor_admitted;
    bool neutral_reference_admitted;
    std::vector<std::string> component_ids;
    std::vector<const char*> component_id_pointers;
    std::vector<int32_t> component_charges;
    std::array<epcsaft_native_capability_coordinate_v1, 4>
        neutral_binary_kij_coordinates;
    std::array<epcsaft_native_capability_coordinate_v1, 4>
        neutral_binary_lij_coordinates;
    std::array<epcsaft_native_capability_coordinate_v1, 3>
        pure_segment_count_coordinates;
    std::array<epcsaft_native_capability_coordinate_v1, 3>
        pure_segment_diameter_coordinates;
    std::array<epcsaft_native_capability_coordinate_v1, 3>
        pure_dispersion_energy_coordinates;
    std::array<epcsaft_native_capability_coordinate_v1, 3>
        pure_association_energy_coordinates;
    std::array<epcsaft_native_capability_coordinate_v1, 3>
        pure_association_volume_coordinates;
    std::array<epcsaft_native_capability_coordinate_v1, 1>
        ion_solvation_born_coordinates;
    std::array<epcsaft_native_capability_coordinate_v1, 1>
        ion_solvation_ionic_permittivity_coordinates;
    std::array<epcsaft_native_capability_coordinate_v1, 1>
        ion_solvation_solvent_permittivity_coordinates;
    std::array<epcsaft_native_capability_coordinate_v1, 1>
        aqueous_solvation_factor_coordinates;
    std::array<epcsaft_native_capability_coordinate_v1, 1>
        aqueous_water_cation_kij_coordinates;
    std::array<epcsaft_native_capability_coordinate_v1, 1>
        aqueous_water_anion_kij_coordinates;
    std::array<epcsaft_native_capability_coordinate_v1, 1>
        aqueous_cation_anion_kij_coordinates;
    std::array<epcsaft_native_capability_coordinate_v1, 1>
        ion_fraction_suppression_coordinates;
    std::vector<epcsaft_native_capability_descriptor_v1> capabilities;
    epcsaft_native_sdk_v1 table;

    NativeSdkEndpoint(
        const NativeModel& source_model,
        std::string fingerprint,
        double minimum_temperature,
        double maximum_temperature,
        std::vector<std::string> identifiers,
        double maximum_ion_fraction,
        bool admit_phase_callbacks,
        bool admit_active_born,
        bool admit_aqueous_kij,
        bool admit_ion_fraction_suppression,
        bool admit_ion_solvation_kij,
        double minimum_pressure,
        double maximum_pressure,
        bool admit_neutral_reference,
        bool admit_aqueous_solvation_factor,
        std::string topology_identity
    )
        : model(source_model),
          parameter_fingerprint(std::move(fingerprint)),
          topology_fingerprint(std::move(topology_identity)),
          temperature_min_k(minimum_temperature),
          temperature_max_k(maximum_temperature),
          pressure_min_pa(minimum_pressure),
          pressure_max_pa(maximum_pressure),
          ion_mole_fraction_max(maximum_ion_fraction),
          phase_callbacks_admitted(admit_phase_callbacks),
          active_born_admitted(admit_active_born),
          aqueous_kij_admitted(admit_aqueous_kij),
          ion_fraction_suppression_admitted(admit_ion_fraction_suppression),
          ion_solvation_kij_admitted(admit_ion_solvation_kij),
          aqueous_solvation_factor_admitted(
              admit_aqueous_solvation_factor
          ),
          neutral_reference_admitted(admit_neutral_reference),
          component_ids(std::move(identifiers)),
          component_id_pointers([&] {
              std::vector<const char*> pointers;
              pointers.reserve(component_ids.size());
              for (const auto& identifier : component_ids) {
                  pointers.push_back(identifier.c_str());
              }
              return pointers;
          }()),
          component_charges(
              source_model.charge_numbers().begin(),
              source_model.charge_numbers().end()
          ),
          neutral_binary_kij_coordinates{{
              {
                  sizeof(epcsaft_native_capability_coordinate_v1),
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_AMOUNT_V1,
                  0,
                  -1,
                  -1,
                  "mol",
              },
              {
                  sizeof(epcsaft_native_capability_coordinate_v1),
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_AMOUNT_V1,
                  1,
                  -1,
                  -1,
                  "mol",
              },
              {
                  sizeof(epcsaft_native_capability_coordinate_v1),
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_VOLUME_V1,
                  -1,
                  -1,
                  -1,
                  "m3",
              },
              {
                  sizeof(epcsaft_native_capability_coordinate_v1),
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_BINARY_INTERACTION_KIJ_V1,
                  -1,
                  0,
                  1,
                  "dimensionless",
              },
          }},
          neutral_binary_lij_coordinates{{
              {
                  sizeof(epcsaft_native_capability_coordinate_v1),
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_AMOUNT_V1,
                  0,
                  -1,
                  -1,
                  "mol",
              },
              {
                  sizeof(epcsaft_native_capability_coordinate_v1),
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_AMOUNT_V1,
                  1,
                  -1,
                  -1,
                  "mol",
              },
              {
                  sizeof(epcsaft_native_capability_coordinate_v1),
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_VOLUME_V1,
                  -1,
                  -1,
                  -1,
                  "m3",
              },
              {
                  sizeof(epcsaft_native_capability_coordinate_v1),
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_BINARY_INTERACTION_LIJ_V1,
                  -1,
                  0,
                  1,
                  "dimensionless",
              },
          }},
          pure_segment_count_coordinates{{
              {
                  sizeof(epcsaft_native_capability_coordinate_v1),
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_AMOUNT_V1,
                  0,
                  -1,
                  -1,
                  "mol",
              },
              {
                  sizeof(epcsaft_native_capability_coordinate_v1),
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_VOLUME_V1,
                  -1,
                  -1,
                  -1,
                  "m3",
              },
              {
                  sizeof(epcsaft_native_capability_coordinate_v1),
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_SEGMENT_COUNT_V1,
                  0,
                  -1,
                  -1,
                  "dimensionless",
              },
          }},
          pure_segment_diameter_coordinates(pure_segment_count_coordinates),
          pure_dispersion_energy_coordinates(pure_segment_count_coordinates),
          pure_association_energy_coordinates(pure_segment_count_coordinates),
          pure_association_volume_coordinates(pure_segment_count_coordinates),
          ion_solvation_born_coordinates{{
              {
                  sizeof(epcsaft_native_capability_coordinate_v1),
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_BORN_DIAMETER_V1,
                  1,
                  -1,
                  -1,
                  "angstrom",
              },
          }},
          ion_solvation_ionic_permittivity_coordinates{{
              {
                  sizeof(epcsaft_native_capability_coordinate_v1),
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_IONIC_REGION_RELATIVE_PERMITTIVITY_V1,
                  -1,
                  -1,
                  -1,
                  "dimensionless",
              },
          }},
          ion_solvation_solvent_permittivity_coordinates{{
              {
                  sizeof(epcsaft_native_capability_coordinate_v1),
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_RELATIVE_PERMITTIVITY_V1,
                  0,
                  -1,
                  -1,
                  "dimensionless",
              },
          }},
          aqueous_solvation_factor_coordinates{{
              {
                  sizeof(epcsaft_native_capability_coordinate_v1),
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_SOLVATION_FACTOR_V1,
                  0,
                  -1,
                  -1,
                  "dimensionless",
              },
          }},
          aqueous_water_cation_kij_coordinates{{
              {
                  sizeof(epcsaft_native_capability_coordinate_v1),
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_BINARY_INTERACTION_KIJ_V1,
                  -1,
                  0,
                  1,
                  "dimensionless",
              },
          }},
          aqueous_water_anion_kij_coordinates{{
              {
                  sizeof(epcsaft_native_capability_coordinate_v1),
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_BINARY_INTERACTION_KIJ_V1,
                  -1,
                  0,
                  2,
                  "dimensionless",
              },
          }},
          aqueous_cation_anion_kij_coordinates{{
              {
                  sizeof(epcsaft_native_capability_coordinate_v1),
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_BINARY_INTERACTION_KIJ_V1,
                  -1,
                  1,
                  2,
                  "dimensionless",
              },
          }},
          ion_fraction_suppression_coordinates{{
              {
                  sizeof(epcsaft_native_capability_coordinate_v1),
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_ION_FRACTION_SUPPRESSION_V1,
                  -1,
                  -1,
                  -1,
                  "dimensionless",
              },
          }},
          capabilities([&] {
              std::vector<epcsaft_native_capability_descriptor_v1> result;
              pure_segment_diameter_coordinates[2].kind =
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_SEGMENT_DIAMETER_V1;
              pure_segment_diameter_coordinates[2].unit = "angstrom";
              pure_dispersion_energy_coordinates[2].kind =
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_DISPERSION_ENERGY_OVER_K_V1;
              pure_dispersion_energy_coordinates[2].unit = "kelvin";
              pure_association_energy_coordinates[2].kind =
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_ASSOCIATION_ENERGY_OVER_K_V1;
              pure_association_energy_coordinates[2].component_index = -1;
              pure_association_energy_coordinates[2].unit = "kelvin";
              pure_association_volume_coordinates[2].kind =
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_ASSOCIATION_VOLUME_V1;
              pure_association_volume_coordinates[2].component_index = -1;
              pure_association_volume_coordinates[2].unit = "dimensionless";
              auto descriptor = [&](
                  uint32_t capability,
                  uint32_t parameter_family,
                  uint32_t parameter_identity,
                  uint32_t model_domain,
                  std::size_t state_coordinate_count,
                  const auto& coordinates
              ) {
                  epcsaft_native_capability_descriptor_v1 value{};
                  value.struct_size =
                      sizeof(epcsaft_native_capability_descriptor_v1);
                  value.schema_version =
                      EPCSAFT_NATIVE_CAPABILITY_SCHEMA_VERSION_V1;
                  value.capability = capability;
                  value.parameter_family = parameter_family;
                  value.parameter_identity = parameter_identity;
                  value.observation_contract =
                      EPCSAFT_NATIVE_OBSERVATION_FIXED_COMPOSITION_HELMHOLTZ_PHASE_V1;
                  value.model_domain = model_domain;
                  value.tensor_layout =
                      EPCSAFT_NATIVE_TENSOR_LAYOUT_ROW_MAJOR_V1;
                  value.derivative_order = 2;
                  value.maturity =
                      EPCSAFT_NATIVE_CAPABILITY_DERIVATIVE_READY_V1;
                  value.authority_effect =
                      EPCSAFT_NATIVE_AUTHORITY_EFFECT_NONE_V1;
                  value.unsupported_status =
                      EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1;
                  value.domain_status = EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1;
                  value.state_coordinate_count = state_coordinate_count;
                  value.active_parameter_count = 1;
                  value.coordinate_count = coordinates.size();
                  value.coordinates = coordinates.data();
                  value.component_count = component_ids.size();
                  value.component_ids = component_id_pointers.data();
                  value.temperature_min_k = temperature_min_k;
                  value.temperature_max_k = temperature_max_k;
                  std::snprintf(
                      value.parameter_fingerprint,
                      sizeof(value.parameter_fingerprint),
                      "%s",
                      parameter_fingerprint.c_str()
                  );
                  std::snprintf(
                      value.topology_fingerprint,
                      sizeof(value.topology_fingerprint),
                      "%s",
                      topology_fingerprint.c_str()
                  );
                  std::snprintf(
                      value.helmholtz_basis_id,
                      sizeof(value.helmholtz_basis_id),
                      "%s",
                      EPCSAFT_NATIVE_SDK_V1_HELMHOLTZ_BASIS_ID
                  );
                  return value;
              };
              auto direct_descriptor = [&](
                  uint32_t capability,
                  uint32_t parameter_family,
                  uint32_t observation_contract,
                  uint32_t model_domain,
                  const auto& coordinates
              ) {
                  auto value = descriptor(
                      capability,
                      parameter_family,
                      EPCSAFT_NATIVE_PARAMETER_IDENTITY_COMPONENT_V1,
                      model_domain,
                      0,
                      coordinates
                  );
                  value.observation_contract = observation_contract;
                  value.derivative_order = 1;
                  value.temperature_min_k = 298.15;
                  value.temperature_max_k = 298.15;
                  return value;
              };
              const bool exact_pure =
                  phase_callbacks_admitted
                  && source_model.supports_parameterized_pure_native_sdk_v1()
                  && component_ids.size() == 1
                  && component_charges == std::vector<int32_t>({0})
                  && !topology_fingerprint.empty();
              if (exact_pure) {
                  result.push_back(descriptor(
                      EPCSAFT_NATIVE_CAPABILITY_PURE_SEGMENT_COUNT_HELMHOLTZ_V1,
                      EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_COUNT_V1,
                      EPCSAFT_NATIVE_PARAMETER_IDENTITY_COMPONENT_V1,
                      EPCSAFT_NATIVE_MODEL_DOMAIN_NEUTRAL_NONASSOCIATING_PURE_V1,
                      2,
                      pure_segment_count_coordinates
                  ));
                  result.push_back(descriptor(
                      EPCSAFT_NATIVE_CAPABILITY_PURE_SEGMENT_DIAMETER_HELMHOLTZ_V1,
                      EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_DIAMETER_V1,
                      EPCSAFT_NATIVE_PARAMETER_IDENTITY_COMPONENT_V1,
                      EPCSAFT_NATIVE_MODEL_DOMAIN_NEUTRAL_NONASSOCIATING_PURE_V1,
                      2,
                      pure_segment_diameter_coordinates
                  ));
                  result.push_back(descriptor(
                      EPCSAFT_NATIVE_CAPABILITY_PURE_DISPERSION_ENERGY_HELMHOLTZ_V1,
                      EPCSAFT_NATIVE_PARAMETER_FAMILY_DISPERSION_ENERGY_OVER_K_V1,
                      EPCSAFT_NATIVE_PARAMETER_IDENTITY_COMPONENT_V1,
                      EPCSAFT_NATIVE_MODEL_DOMAIN_NEUTRAL_NONASSOCIATING_PURE_V1,
                      2,
                      pure_dispersion_energy_coordinates
                  ));
              }
              const bool exact_associating_pure =
                  phase_callbacks_admitted
                  && source_model
                      .supports_associating_pure_parameter_native_sdk_v1()
                  && component_ids.size() == 1
                  && component_charges == std::vector<int32_t>({0})
                  && !topology_fingerprint.empty();
              if (exact_associating_pure) {
                  result.push_back(descriptor(
                      EPCSAFT_NATIVE_CAPABILITY_PURE_ASSOCIATION_ENERGY_HELMHOLTZ_V1,
                      EPCSAFT_NATIVE_PARAMETER_FAMILY_ASSOCIATION_ENERGY_OVER_K_V1,
                      EPCSAFT_NATIVE_PARAMETER_IDENTITY_MODEL_V1,
                      EPCSAFT_NATIVE_MODEL_DOMAIN_NEUTRAL_ASSOCIATING_PURE_V1,
                      2,
                      pure_association_energy_coordinates
                  ));
                  result.push_back(descriptor(
                      EPCSAFT_NATIVE_CAPABILITY_PURE_ASSOCIATION_VOLUME_HELMHOLTZ_V1,
                      EPCSAFT_NATIVE_PARAMETER_FAMILY_ASSOCIATION_VOLUME_V1,
                      EPCSAFT_NATIVE_PARAMETER_IDENTITY_MODEL_V1,
                      EPCSAFT_NATIVE_MODEL_DOMAIN_NEUTRAL_ASSOCIATING_PURE_V1,
                      2,
                      pure_association_volume_coordinates
                  ));
              }
              const bool exact_neutral_binary =
                  phase_callbacks_admitted
                  && source_model.supports_mixture_native_sdk_v1()
                  && component_ids.size() == 2
                  && component_charges == std::vector<int32_t>({0, 0})
                  && !topology_fingerprint.empty();
              if (exact_neutral_binary) {
                  result.push_back(descriptor(
                      EPCSAFT_NATIVE_CAPABILITY_NEUTRAL_BINARY_KIJ_HELMHOLTZ_V1,
                      EPCSAFT_NATIVE_PARAMETER_FAMILY_BINARY_INTERACTION_KIJ_V1,
                      EPCSAFT_NATIVE_PARAMETER_IDENTITY_UNORDERED_COMPONENT_PAIR_V1,
                      EPCSAFT_NATIVE_MODEL_DOMAIN_NEUTRAL_NONASSOCIATING_BINARY_V1,
                      3,
                      neutral_binary_kij_coordinates
                  ));
                  result.push_back(descriptor(
                      EPCSAFT_NATIVE_CAPABILITY_NEUTRAL_BINARY_LIJ_HELMHOLTZ_V1,
                      EPCSAFT_NATIVE_PARAMETER_FAMILY_BINARY_INTERACTION_LIJ_V1,
                      EPCSAFT_NATIVE_PARAMETER_IDENTITY_UNORDERED_COMPONENT_PAIR_V1,
                      EPCSAFT_NATIVE_MODEL_DOMAIN_NEUTRAL_NONASSOCIATING_BINARY_V1,
                      3,
                      neutral_binary_lij_coordinates
                  ));
              }
              if (active_born_admitted) {
                  result.push_back(direct_descriptor(
                      EPCSAFT_NATIVE_CAPABILITY_ION_SOLVATION_BORN_V1,
                      EPCSAFT_NATIVE_PARAMETER_FAMILY_BORN_DIAMETER_V1,
                      EPCSAFT_NATIVE_OBSERVATION_ION_SOLVATION_GIBBS_V1,
                      EPCSAFT_NATIVE_MODEL_DOMAIN_AQUEOUS_SINGLE_ION_V1,
                      ion_solvation_born_coordinates
                  ));
                  auto ionic_permittivity = descriptor(
                      EPCSAFT_NATIVE_CAPABILITY_ION_SOLVATION_IONIC_REGION_PERMITTIVITY_V1,
                      EPCSAFT_NATIVE_PARAMETER_FAMILY_IONIC_REGION_RELATIVE_PERMITTIVITY_V1,
                      EPCSAFT_NATIVE_PARAMETER_IDENTITY_MODEL_V1,
                      EPCSAFT_NATIVE_MODEL_DOMAIN_AQUEOUS_SINGLE_ION_V1,
                      0,
                      ion_solvation_ionic_permittivity_coordinates
                  );
                  ionic_permittivity.observation_contract =
                      EPCSAFT_NATIVE_OBSERVATION_ION_SOLVATION_GIBBS_V1;
                  ionic_permittivity.derivative_order = 1;
                  ionic_permittivity.temperature_min_k = 298.15;
                  ionic_permittivity.temperature_max_k = 298.15;
                  result.push_back(ionic_permittivity);
                  auto solvent_permittivity = descriptor(
                      EPCSAFT_NATIVE_CAPABILITY_ION_SOLVATION_SOLVENT_PERMITTIVITY_V1,
                      EPCSAFT_NATIVE_PARAMETER_FAMILY_RELATIVE_PERMITTIVITY_V1,
                      EPCSAFT_NATIVE_PARAMETER_IDENTITY_COMPONENT_V1,
                      EPCSAFT_NATIVE_MODEL_DOMAIN_AQUEOUS_SINGLE_ION_V1,
                      0,
                      ion_solvation_solvent_permittivity_coordinates
                  );
                  solvent_permittivity.observation_contract =
                      EPCSAFT_NATIVE_OBSERVATION_ION_SOLVATION_GIBBS_V1;
                  solvent_permittivity.derivative_order = 1;
                  solvent_permittivity.temperature_min_k = 298.15;
                  solvent_permittivity.temperature_max_k = 298.15;
                  result.push_back(solvent_permittivity);
              }
              if (aqueous_solvation_factor_admitted) {
                  result.push_back(direct_descriptor(
                      EPCSAFT_NATIVE_CAPABILITY_AQUEOUS_SOLVATION_FACTOR_MIAC_V1,
                      EPCSAFT_NATIVE_PARAMETER_FAMILY_SOLVATION_FACTOR_V1,
                      EPCSAFT_NATIVE_OBSERVATION_AQUEOUS_MEAN_IONIC_ACTIVITY_V1,
                      EPCSAFT_NATIVE_MODEL_DOMAIN_AQUEOUS_NABR_V1,
                      aqueous_solvation_factor_coordinates
                  ));
              }
              if (aqueous_kij_admitted) {
                  const auto aqueous_kij_descriptor = [&](
                      uint32_t capability,
                      const auto& coordinates
                  ) {
                      auto value = descriptor(
                          capability,
                          EPCSAFT_NATIVE_PARAMETER_FAMILY_BINARY_INTERACTION_KIJ_V1,
                          EPCSAFT_NATIVE_PARAMETER_IDENTITY_UNORDERED_COMPONENT_PAIR_V1,
                          EPCSAFT_NATIVE_MODEL_DOMAIN_AQUEOUS_NABR_V1,
                          0,
                          coordinates
                      );
                      value.observation_contract =
                          EPCSAFT_NATIVE_OBSERVATION_AQUEOUS_MEAN_IONIC_ACTIVITY_V1;
                      value.derivative_order = 1;
                      value.temperature_min_k = 298.15;
                      value.temperature_max_k = 298.15;
                      return value;
                  };
                  result.push_back(aqueous_kij_descriptor(
                      EPCSAFT_NATIVE_CAPABILITY_AQUEOUS_WATER_CATION_KIJ_MIAC_V1,
                      aqueous_water_cation_kij_coordinates
                  ));
                  result.push_back(aqueous_kij_descriptor(
                      EPCSAFT_NATIVE_CAPABILITY_AQUEOUS_WATER_ANION_KIJ_MIAC_V1,
                      aqueous_water_anion_kij_coordinates
                  ));
                  result.push_back(aqueous_kij_descriptor(
                      EPCSAFT_NATIVE_CAPABILITY_AQUEOUS_CATION_ANION_KIJ_MIAC_V1,
                      aqueous_cation_anion_kij_coordinates
                  ));
              }
              if (ion_fraction_suppression_admitted) {
                  auto dielectric = descriptor(
                      EPCSAFT_NATIVE_CAPABILITY_ION_FRACTION_SUPPRESSION_V1,
                      EPCSAFT_NATIVE_PARAMETER_FAMILY_ION_FRACTION_SUPPRESSION_V1,
                      EPCSAFT_NATIVE_PARAMETER_IDENTITY_MODEL_V1,
                      EPCSAFT_NATIVE_MODEL_DOMAIN_ION_FRACTION_SUPPRESSION_V1,
                      0,
                      ion_fraction_suppression_coordinates
                  );
                  dielectric.observation_contract =
                      EPCSAFT_NATIVE_OBSERVATION_RELATIVE_PERMITTIVITY_RATIO_V1;
                  dielectric.derivative_order = 1;
                  dielectric.temperature_min_k = 298.15;
                  dielectric.temperature_max_k = 298.15;
                  result.push_back(dielectric);
              }
              if (ion_solvation_kij_admitted) {
                  const auto ion_solvation_kij_descriptor = [&](
                      uint32_t capability,
                      const auto& coordinates
                  ) {
                      auto value = descriptor(
                          capability,
                          EPCSAFT_NATIVE_PARAMETER_FAMILY_BINARY_INTERACTION_KIJ_V1,
                          EPCSAFT_NATIVE_PARAMETER_IDENTITY_UNORDERED_COMPONENT_PAIR_V1,
                          EPCSAFT_NATIVE_MODEL_DOMAIN_ORGANIC_SINGLE_ION_SOLVATION_V1,
                          0,
                          coordinates
                      );
                      value.observation_contract =
                          EPCSAFT_NATIVE_OBSERVATION_ION_SOLVATION_GIBBS_V1;
                      value.derivative_order = 1;
                      value.temperature_min_k = 298.15;
                      value.temperature_max_k = 298.15;
                      return value;
                  };
                  result.push_back(ion_solvation_kij_descriptor(
                      EPCSAFT_NATIVE_CAPABILITY_ION_SOLVATION_SOLVENT_CATION_KIJ_V1,
                      aqueous_water_cation_kij_coordinates
                  ));
                  result.push_back(ion_solvation_kij_descriptor(
                      EPCSAFT_NATIVE_CAPABILITY_ION_SOLVATION_SOLVENT_ANION_KIJ_V1,
                      aqueous_water_anion_kij_coordinates
                  ));
                  result.push_back(ion_solvation_kij_descriptor(
                      EPCSAFT_NATIVE_CAPABILITY_ION_SOLVATION_CATION_ANION_KIJ_V1,
                      aqueous_cation_anion_kij_coordinates
                  ));
              }
              return result;
          }()),
          table{
              EPCSAFT_NATIVE_SDK_V1_ABI_VERSION,
              sizeof(epcsaft_native_sdk_v1),
              sizeof(epcsaft_phase_block_result_v1),
              this,
              &NativeSdkEndpoint::evaluate,
              sizeof(epcsaft_parameterized_phase_block_result_v1),
              &NativeSdkEndpoint::evaluate_parameters,
              source_model.component_count(),
              sizeof(epcsaft_mixture_phase_block_result_v1),
              &NativeSdkEndpoint::evaluate_mixture,
              &NativeSdkEndpoint::evaluate_mixture_kij,
              component_id_pointers.data(),
              component_charges.data(),
              &NativeSdkEndpoint::evaluate_electrolyte,
              &NativeSdkEndpoint::evaluate_molar_volume_bounds,
              &NativeSdkEndpoint::evaluate_packing_fraction,
              temperature_min_k,
              temperature_max_k,
              ion_mole_fraction_max,
              sizeof(epcsaft_ion_solvation_born_result_v1),
              &NativeSdkEndpoint::evaluate_ion_solvation_born,
              sizeof(epcsaft_aqueous_miac_kij_result_v1),
              &NativeSdkEndpoint::evaluate_aqueous_miac_kij,
              neutral_reference_admitted ? component_ids.size() - 1u : 0u,
              sizeof(epcsaft_neutral_reference_result_v1),
              &NativeSdkEndpoint::evaluate_neutral_reference,
              sizeof(epcsaft_aqueous_miac_solvation_factor_result_v1),
              &NativeSdkEndpoint::evaluate_aqueous_miac_solvation_factor,
              &NativeSdkEndpoint::evaluate_aqueous_miac_kij_batch,
              &NativeSdkEndpoint::evaluate_aqueous_miac_solvation_factor_batch,
              sizeof(epcsaft_native_evaluation_budget_v1),
              &NativeSdkEndpoint::evaluate_aqueous_miac_kij_batch_bounded,
              &NativeSdkEndpoint::evaluate_aqueous_miac_solvation_factor_batch_bounded,
              capabilities.size(),
              capabilities.empty() ? nullptr : capabilities.data(),
              &NativeSdkEndpoint::evaluate_mixture_lij,
              &NativeSdkEndpoint::evaluate_pure_parameter,
              sizeof(epcsaft_ion_fraction_suppression_result_v1),
              &NativeSdkEndpoint::evaluate_ion_fraction_suppression,
              sizeof(epcsaft_ion_solvation_kij_result_v1),
              &NativeSdkEndpoint::evaluate_ion_solvation_kij,
              sizeof(epcsaft_ion_solvation_ionic_permittivity_result_v1),
              &NativeSdkEndpoint::evaluate_ion_solvation_ionic_permittivity,
              sizeof(epcsaft_ion_solvation_solvent_permittivity_result_v1),
              &NativeSdkEndpoint::evaluate_ion_solvation_solvent_permittivity,
              sizeof(epcsaft_associating_parameterized_phase_block_result_v1),
              &NativeSdkEndpoint::evaluate_associating_parameters,
          } {}

    // EqExclude: transport and buffer mechanics
    static int evaluate_ion_solvation_ionic_permittivity(
        const void* context,
        const char* expected_fingerprint,
        double temperature_k,
        double pressure_pa,
        double ionic_region_relative_permittivity,
        epcsaft_ion_solvation_ionic_permittivity_result_v1* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if (result->struct_size
            != sizeof(epcsaft_ion_solvation_ionic_permittivity_result_v1)) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        const uint32_t result_size = result->struct_size;
        std::memset(result, 0, sizeof(*result));
        result->struct_size = result_size;
        result->solvation_gibbs_j_per_mol =
            std::numeric_limits<double>::quiet_NaN();
        result->derivative_j_per_mol =
            std::numeric_limits<double>::quiet_NaN();
        result->reference_molality_mol_per_kg =
            std::numeric_limits<double>::quiet_NaN();
        result->reference_convergence_error =
            std::numeric_limits<double>::quiet_NaN();
        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr || expected_fingerprint == nullptr
            || !std::isfinite(temperature_k) || !std::isfinite(pressure_pa)
            || !std::isfinite(ionic_region_relative_permittivity)) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "invalid active ionic-region permittivity arguments"
            );
        }
        const auto& endpoint =
            *static_cast<const NativeSdkEndpoint*>(context);
        if (!endpoint.active_born_admitted
            || std::strcmp(
                endpoint.parameter_fingerprint.c_str(),
                expected_fingerprint
            ) != 0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1,
                "unsupported active ionic-region permittivity model"
            );
        }
        if (temperature_k != 298.15 || pressure_pa != 100000.0
            || ionic_region_relative_permittivity <= 1.0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1,
                "active ionic-region permittivity state is outside the admitted domain"
            );
        }
        try {
            const auto value =
                endpoint.model.evaluate_ion_solvation_ionic_permittivity(
                    temperature_k,
                    pressure_pa,
                    ionic_region_relative_permittivity
                );
            result->solvation_gibbs_j_per_mol =
                value.solvation_gibbs_j_per_mol;
            result->derivative_j_per_mol = value.derivative_j_per_mol;
            result->reference_molality_mol_per_kg =
                value.reference_molality_mol_per_kg;
            result->reference_convergence_error =
                value.reference_convergence_error;
            std::snprintf(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                "%s",
                endpoint.parameter_fingerprint.c_str()
            );
            result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
            return result->status;
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::logic_error& error) {
            return fail(EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1, error.what());
        } catch (const std::exception& error) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1, error.what());
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1,
                "unknown active ionic-region permittivity evaluation failure"
            );
        }
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_ion_solvation_solvent_permittivity(
        const void* context,
        const char* expected_fingerprint,
        double temperature_k,
        double pressure_pa,
        double solvent_relative_permittivity,
        epcsaft_ion_solvation_solvent_permittivity_result_v1* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if (result->struct_size
            != sizeof(epcsaft_ion_solvation_solvent_permittivity_result_v1)) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        const uint32_t result_size = result->struct_size;
        std::memset(result, 0, sizeof(*result));
        result->struct_size = result_size;
        result->solvation_gibbs_j_per_mol =
            std::numeric_limits<double>::quiet_NaN();
        result->derivative_j_per_mol =
            std::numeric_limits<double>::quiet_NaN();
        result->reference_molality_mol_per_kg =
            std::numeric_limits<double>::quiet_NaN();
        result->reference_convergence_error =
            std::numeric_limits<double>::quiet_NaN();
        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr || expected_fingerprint == nullptr
            || !std::isfinite(temperature_k) || !std::isfinite(pressure_pa)
            || !std::isfinite(solvent_relative_permittivity)) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "invalid active solvent permittivity arguments"
            );
        }
        const auto& endpoint =
            *static_cast<const NativeSdkEndpoint*>(context);
        if (!endpoint.active_born_admitted
            || std::strcmp(
                endpoint.parameter_fingerprint.c_str(),
                expected_fingerprint
            ) != 0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1,
                "unsupported active solvent permittivity model"
            );
        }
        if (temperature_k != 298.15 || pressure_pa != 100000.0
            || solvent_relative_permittivity <= 1.0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1,
                "active solvent permittivity state is outside the admitted domain"
            );
        }
        try {
            const auto value =
                endpoint.model.evaluate_ion_solvation_solvent_permittivity(
                    temperature_k,
                    pressure_pa,
                    solvent_relative_permittivity
                );
            result->solvation_gibbs_j_per_mol =
                value.solvation_gibbs_j_per_mol;
            result->derivative_j_per_mol = value.derivative_j_per_mol;
            result->reference_molality_mol_per_kg =
                value.reference_molality_mol_per_kg;
            result->reference_convergence_error =
                value.reference_convergence_error;
            std::snprintf(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                "%s",
                endpoint.parameter_fingerprint.c_str()
            );
            result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
            return result->status;
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::logic_error& error) {
            return fail(EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1, error.what());
        } catch (const std::exception& error) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1, error.what());
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1,
                "unknown active solvent permittivity evaluation failure"
            );
        }
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_ion_solvation_kij(
        const void* context,
        const char* expected_fingerprint,
        double temperature_k,
        double pressure_pa,
        size_t active_ion_component_index,
        const double* k_ij,
        size_t k_ij_count,
        epcsaft_ion_solvation_kij_result_v1* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if (result->struct_size
            != sizeof(epcsaft_ion_solvation_kij_result_v1)) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        const uint32_t result_size = result->struct_size;
        std::memset(result, 0, sizeof(*result));
        result->struct_size = result_size;
        result->solvation_gibbs_j_per_mol =
            std::numeric_limits<double>::quiet_NaN();
        std::fill(
            std::begin(result->derivative),
            std::end(result->derivative),
            std::numeric_limits<double>::quiet_NaN()
        );
        result->reference_molality_mol_per_kg =
            std::numeric_limits<double>::quiet_NaN();
        result->reference_convergence_error =
            std::numeric_limits<double>::quiet_NaN();
        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr || expected_fingerprint == nullptr
            || k_ij == nullptr || k_ij_count != 3
            || !std::isfinite(temperature_k)
            || !std::isfinite(pressure_pa)
            || !std::all_of(k_ij, k_ij + k_ij_count, [](double value) {
                return std::isfinite(value);
            })) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "invalid active ion-solvation k_ij arguments"
            );
        }
        const auto& endpoint =
            *static_cast<const NativeSdkEndpoint*>(context);
        if (!endpoint.ion_solvation_kij_admitted
            || std::strcmp(
                endpoint.parameter_fingerprint.c_str(),
                expected_fingerprint
            ) != 0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1,
                "unsupported active ion-solvation k_ij model"
            );
        }
        if (temperature_k != 298.15 || pressure_pa != 100000.0
            || (active_ion_component_index != 1
                && active_ion_component_index != 2)) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1,
                "active ion-solvation k_ij state is outside the admitted domain"
            );
        }
        try {
            const std::array<double, 3> coordinates{
                k_ij[0], k_ij[1], k_ij[2]
            };
            const auto value = endpoint.model.evaluate_ion_solvation_kij(
                temperature_k,
                pressure_pa,
                active_ion_component_index,
                coordinates
            );
            result->solvation_gibbs_j_per_mol =
                value.solvation_gibbs_j_per_mol;
            std::copy(
                value.derivative.begin(),
                value.derivative.end(),
                result->derivative
            );
            result->reference_molality_mol_per_kg =
                value.reference_molality_mol_per_kg;
            result->reference_convergence_error =
                value.reference_convergence_error;
            std::snprintf(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                "%s",
                endpoint.parameter_fingerprint.c_str()
            );
            result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
            return result->status;
        } catch (const detail::DensityDomainError& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::logic_error& error) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1, error.what()
            );
        } catch (const std::exception& error) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1, error.what()
            );
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1,
                "unknown active ion-solvation k_ij failure"
            );
        }
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_ion_fraction_suppression(
        const void* context,
        const char* expected_fingerprint,
        double temperature_k,
        double pressure_pa,
        double total_ion_mole_fraction,
        double suppression_coefficient,
        epcsaft_ion_fraction_suppression_result_v1* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if (result->struct_size
            != sizeof(epcsaft_ion_fraction_suppression_result_v1)) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        const double nan = std::numeric_limits<double>::quiet_NaN();
        result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
        result->bulk_relative_permittivity = nan;
        result->derivative = nan;
        result->salt_free_relative_permittivity = nan;
        result->parameter_fingerprint[0] = '\0';
        result->error[0] = '\0';
        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr || expected_fingerprint == nullptr
            || !std::isfinite(temperature_k)
            || !std::isfinite(pressure_pa)
            || temperature_k <= 0.0 || pressure_pa <= 0.0
            || !std::isfinite(total_ion_mole_fraction)
            || !std::isfinite(suppression_coefficient)) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "invalid active ion-fraction suppression arguments"
            );
        }
        const auto& endpoint =
            *static_cast<const NativeSdkEndpoint*>(context);
        if (!endpoint.ion_fraction_suppression_admitted
            || std::strcmp(
                expected_fingerprint,
                endpoint.parameter_fingerprint.c_str()
            ) != 0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1,
                "unsupported active ion-fraction suppression model"
            );
        }
        if (temperature_k != 298.15 || pressure_pa != 100000.0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1,
                "state is outside the active ion-fraction suppression domain"
            );
        }
        try {
            const auto value = endpoint.model.evaluate_ion_fraction_suppression(
                total_ion_mole_fraction,
                suppression_coefficient
            );
            result->bulk_relative_permittivity =
                value.bulk_relative_permittivity;
            result->derivative = value.derivative;
            result->salt_free_relative_permittivity =
                value.salt_free_relative_permittivity;
            std::snprintf(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                "%s",
                endpoint.parameter_fingerprint.c_str()
            );
            return result->status;
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::logic_error& error) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1, error.what()
            );
        } catch (const std::exception& error) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1, error.what()
            );
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1,
                "unknown active ion-fraction suppression failure"
            );
        }
    }

    // EqID: neutral_reference_stoichiometric_basis
    static int evaluate_neutral_reference(
        const void* context,
        const char* expected_fingerprint,
        double temperature_k,
        double pressure_pa,
        epcsaft_neutral_reference_result_v1* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if (result->struct_size != sizeof(epcsaft_neutral_reference_result_v1)) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        const double nan = std::numeric_limits<double>::quiet_NaN();
        result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
        result->derivative_availability =
            EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_NONE_V1;
        result->temperature_k = nan;
        result->pressure_pa = nan;
        result->solvent_molar_mass_kg_per_mol = nan;
        result->reference_amount_mol = nan;
        result->reference_number_density_mol_per_m3 = nan;
        result->reference_molality_mol_per_kg = nan;
        result->reference_convergence_error = nan;
        result->parameter_fingerprint[0] = '\0';
        result->helmholtz_basis_id[0] = '\0';
        result->error[0] = '\0';
        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr || expected_fingerprint == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "model context and expected fingerprint are required"
            );
        }
        if (!std::isfinite(temperature_k) || temperature_k <= 0.0
            || !std::isfinite(pressure_pa) || pressure_pa <= 0.0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "temperature and pressure must be positive and finite"
            );
        }
        const auto& endpoint = *static_cast<const NativeSdkEndpoint*>(context);
        if (!endpoint.neutral_reference_admitted
            || std::strcmp(
                endpoint.parameter_fingerprint.c_str(), expected_fingerprint
            ) != 0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1,
                "unsupported neutral standard-reference model identity"
            );
        }
        if (!std::isfinite(endpoint.temperature_min_k)
            || !std::isfinite(endpoint.temperature_max_k)
            || !std::isfinite(endpoint.pressure_min_pa)
            || !std::isfinite(endpoint.pressure_max_pa)
            || temperature_k < endpoint.temperature_min_k
            || temperature_k > endpoint.temperature_max_k
            || pressure_pa < endpoint.pressure_min_pa
            || pressure_pa > endpoint.pressure_max_pa) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1,
                "standard-reference state is outside the parameter source domain"
            );
        }
        const std::size_t component_count = endpoint.component_ids.size();
        const std::size_t row_count = component_count - 1;
        if (result->component_count != component_count
            || result->neutral_basis_row_count != row_count
            || result->neutral_basis_capacity != row_count * component_count
            || result->contraction_capacity != row_count
            || result->reference_composition_capacity != component_count
            || result->neutral_basis == nullptr
            || result->log_fugacity_contractions == nullptr
            || result->reference_composition == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "neutral reference buffers do not match the immutable dimensions"
            );
        }
        try {
            std::size_t solvent = component_count;
            std::vector<std::size_t> cations;
            std::vector<std::size_t> anions;
            for (std::size_t component = 0; component < component_count; ++component) {
                if (endpoint.component_ids[component] == "water"
                    && endpoint.component_charges[component] == 0) {
                    solvent = component;
                }
                if (endpoint.component_charges[component] > 0) {
                    cations.push_back(component);
                } else if (endpoint.component_charges[component] < 0) {
                    anions.push_back(component);
                }
            }
            if (solvent == component_count || cations.empty() || anions.empty()) {
                return fail(
                    EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1,
                    "neutral reference component topology is incomplete"
                );
            }
            std::vector<double> reference_composition(component_count, 0.0);
            reference_composition[solvent] = 1.0;
            const auto pure_state = endpoint.model.evaluate_pressure(
                temperature_k, pressure_pa, reference_composition, "liquid"
            );
            if (!pure_state.log_fugacity_coefficient.has_value()) {
                throw std::runtime_error("pure-solvent reference fugacity is unavailable");
            }
            std::vector<double> basis(row_count * component_count, 0.0);
            std::vector<double> contractions(row_count, 0.0);
            std::size_t row = 0;
            for (std::size_t component = 0; component < component_count; ++component) {
                if (endpoint.component_charges[component] == 0) {
                    basis[row * component_count + component] = 1.0;
                    contractions[row] = (*pure_state.log_fugacity_coefficient)[component];
                    ++row;
                }
            }
            double reference_molality = std::numeric_limits<double>::quiet_NaN();
            double reference_convergence_error = 0.0;
            auto append_pair = [&](std::size_t cation, std::size_t anion) {
                const int divisor = std::gcd(
                    std::abs(endpoint.component_charges[cation]),
                    std::abs(endpoint.component_charges[anion])
                );
                const std::array<int, 2> stoichiometry = {
                    std::abs(endpoint.component_charges[anion]) / divisor,
                    std::abs(endpoint.component_charges[cation]) / divisor,
                };
                const auto sequence = endpoint.model.electrolyte_reference_sequence(
                    temperature_k,
                    pressure_pa,
                    reference_composition,
                    std::array<std::size_t, 2>{cation, anion},
                    stoichiometry,
                    "liquid"
                );
                basis[row * component_count + cation] = stoichiometry[0];
                basis[row * component_count + anion] = stoichiometry[1];
                contractions[row] = stoichiometry[0] * sequence.log_fugacity[0]
                    + stoichiometry[1] * sequence.log_fugacity[1];
                if (std::isnan(reference_molality)) {
                    reference_molality = sequence.reference_molality;
                } else if (sequence.reference_molality != reference_molality) {
                    throw std::runtime_error("neutral reference molality changed between rows");
                }
                reference_convergence_error = std::max(
                    reference_convergence_error, sequence.convergence_error
                );
                ++row;
            };
            for (std::size_t anion : anions) {
                append_pair(cations.front(), anion);
            }
            for (std::size_t index = 1; index < cations.size(); ++index) {
                append_pair(cations[index], anions.front());
            }
            if (row != row_count) {
                throw std::runtime_error("neutral reference basis dimension changed");
            }
            const auto& molar_masses = endpoint.model.molar_masses_kg_per_mol();
            const double solvent_molar_mass = molar_masses[solvent];
            if (!std::all_of(basis.begin(), basis.end(), [](double value) {
                    return std::isfinite(value);
                })
                || !std::all_of(
                    contractions.begin(), contractions.end(), [](double value) {
                        return std::isfinite(value);
                    }
                )
                || !std::isfinite(solvent_molar_mass) || solvent_molar_mass <= 0.0
                || !std::isfinite(reference_molality)
                || reference_molality <= 0.0
                || !std::isfinite(reference_convergence_error)
                || reference_convergence_error < 0.0) {
                throw std::runtime_error("neutral standard-reference result is nonfinite");
            }
            std::copy(basis.begin(), basis.end(), result->neutral_basis);
            std::copy(
                contractions.begin(),
                contractions.end(),
                result->log_fugacity_contractions
            );
            std::copy(
                reference_composition.begin(),
                reference_composition.end(),
                result->reference_composition
            );
            result->temperature_k = temperature_k;
            result->pressure_pa = pressure_pa;
            result->solvent_molar_mass_kg_per_mol = solvent_molar_mass;
            result->reference_amount_mol = 1.0;
            result->reference_number_density_mol_per_m3 = 1.0;
            result->reference_molality_mol_per_kg = reference_molality;
            result->reference_convergence_error = reference_convergence_error;
            std::snprintf(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                "%s",
                endpoint.parameter_fingerprint.c_str()
            );
            std::snprintf(
                result->helmholtz_basis_id,
                sizeof(result->helmholtz_basis_id),
                "%s",
                EPCSAFT_NATIVE_HELMHOLTZ_BASIS_ID_V1
            );
            result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
            return result->status;
        } catch (const detail::DensityDomainError& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::logic_error& error) {
            return fail(EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1, error.what());
        } catch (const std::exception& error) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1, error.what());
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1,
                "unknown neutral standard-reference evaluation failure"
            );
        }
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_aqueous_miac_kij(
        const void* context,
        const char* expected_fingerprint,
        double temperature_k,
        double pressure_pa,
        double molality,
        const double* k_ij,
        size_t k_ij_count,
        epcsaft_aqueous_miac_kij_result_v1* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if (result->struct_size != sizeof(epcsaft_aqueous_miac_kij_result_v1)) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        const uint32_t result_size = result->struct_size;
        std::memset(result, 0, sizeof(*result));
        result->struct_size = result_size;
        result->log_mean_ionic_activity_coefficient_molality =
            std::numeric_limits<double>::quiet_NaN();
        std::fill(
            std::begin(result->derivative),
            std::end(result->derivative),
            std::numeric_limits<double>::quiet_NaN()
        );
        result->reference_molality_mol_per_kg =
            std::numeric_limits<double>::quiet_NaN();
        result->reference_convergence_error =
            std::numeric_limits<double>::quiet_NaN();
        result->reference_derivative_convergence_error =
            std::numeric_limits<double>::quiet_NaN();
        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr || expected_fingerprint == nullptr || k_ij == nullptr) {
            return fail(EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1, "required pointer is null");
        }
        if (k_ij_count != EPCSAFT_AQUEOUS_ACTIVE_KIJ_COORDINATE_COUNT_V1) {
            return fail(EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1, "active k_ij count must be three");
        }
        if (!std::isfinite(temperature_k) || !std::isfinite(pressure_pa)
            || !std::isfinite(molality)
            || !std::all_of(k_ij, k_ij + k_ij_count, [](double value) {
                return std::isfinite(value);
            })) {
            return fail(EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1, "state and active k_ij must be finite");
        }
        const auto& endpoint = *static_cast<const NativeSdkEndpoint*>(context);
        if (!endpoint.aqueous_kij_admitted
            || std::strcmp(
                endpoint.parameter_fingerprint.c_str(), expected_fingerprint
            ) != 0) {
            return fail(EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1, "unsupported aqueous active k_ij model");
        }
        if (temperature_k != 298.15 || pressure_pa != 100000.0
            || molality < 0.001 || molality > 6.0) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, "aqueous active k_ij state is outside the admitted domain");
        }
        try {
            const std::array<double, 3> coordinates = {k_ij[0], k_ij[1], k_ij[2]};
            const auto value = endpoint.model.evaluate_aqueous_miac_kij(
                temperature_k, pressure_pa, molality, coordinates
            );
            result->log_mean_ionic_activity_coefficient_molality =
                value.log_mean_ionic_activity_coefficient_molality;
            std::copy(value.derivative.begin(), value.derivative.end(), result->derivative);
            result->reference_molality_mol_per_kg = value.reference_molality_mol_per_kg;
            result->reference_convergence_error = value.reference_convergence_error;
            result->reference_derivative_convergence_error =
                value.reference_derivative_convergence_error;
            std::snprintf(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                "%s",
                endpoint.parameter_fingerprint.c_str()
            );
            result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
            return result->status;
        } catch (const detail::DensityDomainError& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::logic_error& error) {
            return fail(EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1, error.what());
        } catch (const std::exception& error) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1, error.what());
        } catch (...) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1, "unknown aqueous active k_ij evaluation failure");
        }
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_aqueous_miac_solvation_factor(
        const void* context,
        const char* expected_fingerprint,
        double temperature_k,
        double pressure_pa,
        double molality,
        double solvation_factor,
        epcsaft_aqueous_miac_solvation_factor_result_v1* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if (result->struct_size
            != sizeof(epcsaft_aqueous_miac_solvation_factor_result_v1)) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        const uint32_t result_size = result->struct_size;
        std::memset(result, 0, sizeof(*result));
        result->struct_size = result_size;
        result->log_mean_ionic_activity_coefficient_molality =
            std::numeric_limits<double>::quiet_NaN();
        result->derivative = std::numeric_limits<double>::quiet_NaN();
        result->reference_molality_mol_per_kg =
            std::numeric_limits<double>::quiet_NaN();
        result->reference_convergence_error =
            std::numeric_limits<double>::quiet_NaN();
        result->reference_derivative_convergence_error =
            std::numeric_limits<double>::quiet_NaN();
        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr || expected_fingerprint == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "model context and expected fingerprint are required"
            );
        }
        if (!std::isfinite(temperature_k) || !std::isfinite(pressure_pa)
            || !std::isfinite(molality) || !std::isfinite(solvation_factor)) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "state and active solvation factor must be finite"
            );
        }
        const auto& endpoint = *static_cast<const NativeSdkEndpoint*>(context);
        if (!endpoint.aqueous_solvation_factor_admitted
            || std::strcmp(
                endpoint.parameter_fingerprint.c_str(), expected_fingerprint
            ) != 0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1,
                "unsupported aqueous active solvation-factor model"
            );
        }
        if (temperature_k != 298.15 || pressure_pa != 100000.0
            || molality < 0.001 || molality > 6.0
            || solvation_factor <= 0.0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1,
                "aqueous active solvation-factor state is outside the admitted domain"
            );
        }
        try {
            const auto value =
                endpoint.model.evaluate_aqueous_miac_solvation_factor(
                    temperature_k,
                    pressure_pa,
                    molality,
                    solvation_factor
                );
            result->log_mean_ionic_activity_coefficient_molality =
                value.log_mean_ionic_activity_coefficient_molality;
            result->derivative = value.derivative;
            result->reference_molality_mol_per_kg =
                value.reference_molality_mol_per_kg;
            result->reference_convergence_error =
                value.reference_convergence_error;
            result->reference_derivative_convergence_error =
                value.reference_derivative_convergence_error;
            std::snprintf(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                "%s",
                endpoint.parameter_fingerprint.c_str()
            );
            result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
            return result->status;
        } catch (const detail::DensityDomainError& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::logic_error& error) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1, error.what()
            );
        } catch (const std::exception& error) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1, error.what());
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1,
                "unknown aqueous active solvation-factor evaluation failure"
            );
        }
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_aqueous_miac_kij_batch_impl(
        const void* context,
        const char* expected_fingerprint,
        double temperature_k,
        double pressure_pa,
        const double* molalities,
        size_t molality_count,
        const double* k_ij,
        size_t k_ij_count,
        epcsaft_aqueous_miac_kij_result_v1* results,
        size_t result_count,
        const epcsaft_native_evaluation_budget_v1* budget
    ) noexcept {
        if (results == nullptr || result_count == 0) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        for (std::size_t index = 0; index < result_count; ++index) {
            if (results[index].struct_size
                != sizeof(epcsaft_aqueous_miac_kij_result_v1)) {
                return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
            }
        }
        const double nan = std::numeric_limits<double>::quiet_NaN();
        for (std::size_t index = 0; index < result_count; ++index) {
            const uint32_t result_size = results[index].struct_size;
            std::memset(&results[index], 0, sizeof(results[index]));
            results[index].struct_size = result_size;
            results[index].log_mean_ionic_activity_coefficient_molality = nan;
            std::fill(
                std::begin(results[index].derivative),
                std::end(results[index].derivative),
                nan
            );
            results[index].reference_molality_mol_per_kg = nan;
            results[index].reference_convergence_error = nan;
            results[index].reference_derivative_convergence_error = nan;
        }
        auto fail = [&](int status, const char* message) {
            for (std::size_t index = 0; index < result_count; ++index) {
                results[index].status = status;
                std::snprintf(
                    results[index].error,
                    sizeof(results[index].error),
                    "%s",
                    message
                );
            }
            return status;
        };
        if (context == nullptr || expected_fingerprint == nullptr
            || molalities == nullptr || k_ij == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "required pointer is null"
            );
        }
        if (molality_count != result_count
            || k_ij_count != EPCSAFT_AQUEOUS_ACTIVE_KIJ_COORDINATE_COUNT_V1) {
            return fail(
                EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1,
                "batch dimensions do not match"
            );
        }
        if (!std::isfinite(temperature_k) || !std::isfinite(pressure_pa)
            || !std::all_of(
                molalities,
                molalities + molality_count,
                [](double value) { return std::isfinite(value); }
            )
            || !std::all_of(k_ij, k_ij + k_ij_count, [](double value) {
                return std::isfinite(value);
            })) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "state and active k_ij must be finite"
            );
        }
        const auto& endpoint = *static_cast<const NativeSdkEndpoint*>(context);
        if (!endpoint.aqueous_kij_admitted
            || std::strcmp(
                endpoint.parameter_fingerprint.c_str(), expected_fingerprint
            ) != 0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1,
                "unsupported aqueous active k_ij model"
            );
        }
        if (temperature_k != 298.15 || pressure_pa != 100000.0
            || std::any_of(
                molalities,
                molalities + molality_count,
                [](double value) { return value < 0.001 || value > 6.0; }
            )) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1,
                "aqueous active k_ij state is outside the admitted domain"
            );
        }
        std::optional<NativeEvaluationBudget> bounded_budget;
        if (budget != nullptr) {
            if (budget->struct_size
                != sizeof(epcsaft_native_evaluation_budget_v1)) {
                return fail(
                    EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1,
                    "evaluation budget size does not match"
                );
            }
            bounded_budget.emplace(budget->work_budget, budget->deadline_ns);
        }
        try {
            const std::array<double, 3> coordinates = {
                k_ij[0], k_ij[1], k_ij[2]
            };
            const auto values = endpoint.model.evaluate_aqueous_miac_kij_batch(
                temperature_k,
                pressure_pa,
                std::vector<double>(
                    molalities, molalities + molality_count
                ),
                coordinates,
                bounded_budget.has_value() ? &*bounded_budget : nullptr
            );
            if (values.size() != result_count) {
                throw std::runtime_error(
                    "aqueous active k_ij batch result count mismatch"
                );
            }
            for (std::size_t index = 0; index < result_count; ++index) {
                results[index].log_mean_ionic_activity_coefficient_molality =
                    values[index].log_mean_ionic_activity_coefficient_molality;
                std::copy(
                    values[index].derivative.begin(),
                    values[index].derivative.end(),
                    results[index].derivative
                );
                results[index].reference_molality_mol_per_kg =
                    values[index].reference_molality_mol_per_kg;
                results[index].reference_convergence_error =
                    values[index].reference_convergence_error;
                results[index].reference_derivative_convergence_error =
                    values[index].reference_derivative_convergence_error;
                std::snprintf(
                    results[index].parameter_fingerprint,
                    sizeof(results[index].parameter_fingerprint),
                    "%s",
                    endpoint.parameter_fingerprint.c_str()
                );
                results[index].status = EPCSAFT_NATIVE_STATUS_OK_V1;
            }
            return EPCSAFT_NATIVE_STATUS_OK_V1;
        } catch (const detail::DensityDomainError& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const NativeEvaluationCancelled& error) {
            return fail(EPCSAFT_NATIVE_STATUS_CANCELLED_V1, error.what());
        } catch (const NativeEvaluationDeadlineExceeded& error) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DEADLINE_EXCEEDED_V1, error.what()
            );
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::logic_error& error) {
            return fail(EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1, error.what());
        } catch (const std::exception& error) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1, error.what());
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1,
                "unknown aqueous active k_ij batch evaluation failure"
            );
        }
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_aqueous_miac_kij_batch(
        const void* context,
        const char* expected_fingerprint,
        double temperature_k,
        double pressure_pa,
        const double* molalities,
        size_t molality_count,
        const double* k_ij,
        size_t k_ij_count,
        epcsaft_aqueous_miac_kij_result_v1* results,
        size_t result_count
    ) noexcept {
        return evaluate_aqueous_miac_kij_batch_impl(
            context,
            expected_fingerprint,
            temperature_k,
            pressure_pa,
            molalities,
            molality_count,
            k_ij,
            k_ij_count,
            results,
            result_count,
            nullptr
        );
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_aqueous_miac_kij_batch_bounded(
        const void* context,
        const char* expected_fingerprint,
        double temperature_k,
        double pressure_pa,
        const double* molalities,
        size_t molality_count,
        const double* k_ij,
        size_t k_ij_count,
        epcsaft_aqueous_miac_kij_result_v1* results,
        size_t result_count,
        const epcsaft_native_evaluation_budget_v1* budget
    ) noexcept {
        return evaluate_aqueous_miac_kij_batch_impl(
            context,
            expected_fingerprint,
            temperature_k,
            pressure_pa,
            molalities,
            molality_count,
            k_ij,
            k_ij_count,
            results,
            result_count,
            budget
        );
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_aqueous_miac_solvation_factor_batch_impl(
        const void* context,
        const char* expected_fingerprint,
        double temperature_k,
        double pressure_pa,
        const double* molalities,
        size_t molality_count,
        double solvation_factor,
        epcsaft_aqueous_miac_solvation_factor_result_v1* results,
        size_t result_count,
        const epcsaft_native_evaluation_budget_v1* budget
    ) noexcept {
        if (results == nullptr || result_count == 0) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        for (std::size_t index = 0; index < result_count; ++index) {
            if (results[index].struct_size
                != sizeof(epcsaft_aqueous_miac_solvation_factor_result_v1)) {
                return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
            }
        }
        const double nan = std::numeric_limits<double>::quiet_NaN();
        for (std::size_t index = 0; index < result_count; ++index) {
            const uint32_t result_size = results[index].struct_size;
            std::memset(&results[index], 0, sizeof(results[index]));
            results[index].struct_size = result_size;
            results[index].log_mean_ionic_activity_coefficient_molality = nan;
            results[index].derivative = nan;
            results[index].reference_molality_mol_per_kg = nan;
            results[index].reference_convergence_error = nan;
            results[index].reference_derivative_convergence_error = nan;
        }
        auto fail = [&](int status, const char* message) {
            for (std::size_t index = 0; index < result_count; ++index) {
                results[index].status = status;
                std::snprintf(
                    results[index].error,
                    sizeof(results[index].error),
                    "%s",
                    message
                );
            }
            return status;
        };
        if (context == nullptr || expected_fingerprint == nullptr
            || molalities == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "required pointer is null"
            );
        }
        if (molality_count != result_count) {
            return fail(
                EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1,
                "batch dimensions do not match"
            );
        }
        if (!std::isfinite(temperature_k) || !std::isfinite(pressure_pa)
            || !std::isfinite(solvation_factor)
            || !std::all_of(
                molalities,
                molalities + molality_count,
                [](double value) { return std::isfinite(value); }
            )) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "state and active solvation factor must be finite"
            );
        }
        const auto& endpoint = *static_cast<const NativeSdkEndpoint*>(context);
        if (!endpoint.aqueous_solvation_factor_admitted
            || std::strcmp(
                endpoint.parameter_fingerprint.c_str(), expected_fingerprint
            ) != 0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1,
                "unsupported aqueous active solvation-factor model"
            );
        }
        if (temperature_k != 298.15 || pressure_pa != 100000.0
            || solvation_factor <= 0.0
            || std::any_of(
                molalities,
                molalities + molality_count,
                [](double value) { return value < 0.001 || value > 6.0; }
            )) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1,
                "aqueous active solvation-factor state is outside the admitted domain"
            );
        }
        std::optional<NativeEvaluationBudget> bounded_budget;
        if (budget != nullptr) {
            if (budget->struct_size
                != sizeof(epcsaft_native_evaluation_budget_v1)) {
                return fail(
                    EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1,
                    "evaluation budget size does not match"
                );
            }
            bounded_budget.emplace(budget->work_budget, budget->deadline_ns);
        }
        try {
            const auto values =
                endpoint.model.evaluate_aqueous_miac_solvation_factor_batch(
                    temperature_k,
                    pressure_pa,
                    std::vector<double>(
                        molalities, molalities + molality_count
                    ),
                    solvation_factor,
                    bounded_budget.has_value() ? &*bounded_budget : nullptr
                );
            if (values.size() != result_count) {
                throw std::runtime_error(
                    "aqueous active solvation-factor batch result count mismatch"
                );
            }
            for (std::size_t index = 0; index < result_count; ++index) {
                results[index].log_mean_ionic_activity_coefficient_molality =
                    values[index].log_mean_ionic_activity_coefficient_molality;
                results[index].derivative = values[index].derivative;
                results[index].reference_molality_mol_per_kg =
                    values[index].reference_molality_mol_per_kg;
                results[index].reference_convergence_error =
                    values[index].reference_convergence_error;
                results[index].reference_derivative_convergence_error =
                    values[index].reference_derivative_convergence_error;
                std::snprintf(
                    results[index].parameter_fingerprint,
                    sizeof(results[index].parameter_fingerprint),
                    "%s",
                    endpoint.parameter_fingerprint.c_str()
                );
                results[index].status = EPCSAFT_NATIVE_STATUS_OK_V1;
            }
            return EPCSAFT_NATIVE_STATUS_OK_V1;
        } catch (const detail::DensityDomainError& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const NativeEvaluationCancelled& error) {
            return fail(EPCSAFT_NATIVE_STATUS_CANCELLED_V1, error.what());
        } catch (const NativeEvaluationDeadlineExceeded& error) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DEADLINE_EXCEEDED_V1, error.what()
            );
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::logic_error& error) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1, error.what()
            );
        } catch (const std::exception& error) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1, error.what());
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1,
                "unknown aqueous active solvation-factor batch evaluation failure"
            );
        }
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_aqueous_miac_solvation_factor_batch(
        const void* context,
        const char* expected_fingerprint,
        double temperature_k,
        double pressure_pa,
        const double* molalities,
        size_t molality_count,
        double solvation_factor,
        epcsaft_aqueous_miac_solvation_factor_result_v1* results,
        size_t result_count
    ) noexcept {
        return evaluate_aqueous_miac_solvation_factor_batch_impl(
            context,
            expected_fingerprint,
            temperature_k,
            pressure_pa,
            molalities,
            molality_count,
            solvation_factor,
            results,
            result_count,
            nullptr
        );
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_aqueous_miac_solvation_factor_batch_bounded(
        const void* context,
        const char* expected_fingerprint,
        double temperature_k,
        double pressure_pa,
        const double* molalities,
        size_t molality_count,
        double solvation_factor,
        epcsaft_aqueous_miac_solvation_factor_result_v1* results,
        size_t result_count,
        const epcsaft_native_evaluation_budget_v1* budget
    ) noexcept {
        return evaluate_aqueous_miac_solvation_factor_batch_impl(
            context,
            expected_fingerprint,
            temperature_k,
            pressure_pa,
            molalities,
            molality_count,
            solvation_factor,
            results,
            result_count,
            budget
        );
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_ion_solvation_born(
        const void* context,
        double temperature_k,
        double pressure_pa,
        double born_diameter_angstrom,
        epcsaft_ion_solvation_born_result_v1* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if (result->struct_size
            != sizeof(epcsaft_ion_solvation_born_result_v1)) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        const uint32_t result_size = result->struct_size;
        std::memset(result, 0, sizeof(*result));
        result->struct_size = result_size;
        result->solvation_gibbs_j_per_mol = std::numeric_limits<double>::quiet_NaN();
        result->derivative_j_per_mol_per_angstrom =
            std::numeric_limits<double>::quiet_NaN();
        result->reference_molality_mol_per_kg =
            std::numeric_limits<double>::quiet_NaN();
        result->reference_convergence_error =
            std::numeric_limits<double>::quiet_NaN();
        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "model context is null"
            );
        }
        if (!std::isfinite(temperature_k) || !std::isfinite(pressure_pa)
            || !std::isfinite(born_diameter_angstrom) || temperature_k <= 0.0
            || pressure_pa <= 0.0 || born_diameter_angstrom <= 0.0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "temperature, pressure, and Born diameter must be positive and finite"
            );
        }
        const auto& endpoint = *static_cast<const NativeSdkEndpoint*>(context);
        if (!endpoint.active_born_admitted) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1,
                "model is not admitted for active Born solvation"
            );
        }
        if (temperature_k != 298.15 || pressure_pa != 100000.0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1,
                "active Born reference state requires 298.15 K and 100000 Pa"
            );
        }
        try {
            const auto value = endpoint.model.evaluate_ion_solvation_born(
                temperature_k, pressure_pa, born_diameter_angstrom
            );
            result->solvation_gibbs_j_per_mol = value.solvation_gibbs_j_per_mol;
            result->derivative_j_per_mol_per_angstrom =
                value.derivative_j_per_mol_per_angstrom;
            result->reference_molality_mol_per_kg =
                value.reference_molality_mol_per_kg;
            result->reference_convergence_error =
                value.reference_convergence_error;
            std::snprintf(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                "%s",
                endpoint.parameter_fingerprint.c_str()
            );
            result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
            return result->status;
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::logic_error& error) {
            return fail(EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1, error.what());
        } catch (const std::exception& error) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1, error.what());
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1,
                "unknown active Born evaluation failure"
            );
        }
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate(
        const void* context,
        double temperature_k,
        double amount_mol,
        double volume_m3,
        epcsaft_phase_block_result_v1* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if (result->struct_size != sizeof(epcsaft_phase_block_result_v1)) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        const uint32_t result_size = result->struct_size;
        std::memset(result, 0, sizeof(*result));
        result->struct_size = result_size;

        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "model context is null"
            );
        }
        if (!std::isfinite(temperature_k) || !std::isfinite(amount_mol)
            || !std::isfinite(volume_m3) || temperature_k <= 0.0
            || amount_mol <= 0.0 || volume_m3 <= 0.0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "temperature, amount, and volume must be positive and finite"
            );
        }

        const auto& endpoint = *static_cast<const NativeSdkEndpoint*>(context);
        if (!endpoint.phase_callbacks_admitted) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1,
                "model is not admitted for this phase callback"
            );
        }
        if ((std::isfinite(endpoint.temperature_min_k)
             && temperature_k < endpoint.temperature_min_k)
            || (std::isfinite(endpoint.temperature_max_k)
                && temperature_k > endpoint.temperature_max_k)) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1,
                "temperature is outside the parameter source domain"
            );
        }

        try {
            const auto block = endpoint.model.evaluate_pure_phase(
                temperature_k, amount_mol, volume_m3
            );
            result->helmholtz_over_rt_reference_amount =
                block.helmholtz_over_rt_reference_amount;
            std::copy(
                block.gradient.begin(), block.gradient.end(), result->gradient
            );
            std::copy(
                block.hessian.begin(), block.hessian.end(), result->hessian
            );
            std::copy(block.third.begin(), block.third.end(), result->third);
            result->pressure_pa = block.pressure_pa;
            result->chemical_potential_over_rt = block.chemical_potential_over_rt;
            std::snprintf(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                "%s",
                endpoint.parameter_fingerprint.c_str()
            );
            result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
            return result->status;
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::logic_error& error) {
            return fail(EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1, error.what());
        } catch (const std::exception& error) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1, error.what());
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1,
                "unknown provider evaluation failure"
            );
        }
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_parameters(
        const void* context,
        double temperature_k,
        double amount_mol,
        double volume_m3,
        double segment_count,
        double segment_diameter_angstrom,
        double dispersion_energy_over_k,
        epcsaft_parameterized_phase_block_result_v1* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if (result->struct_size
            != sizeof(epcsaft_parameterized_phase_block_result_v1)) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        const uint32_t result_size = result->struct_size;
        std::memset(result, 0, sizeof(*result));
        result->struct_size = result_size;

        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "model context is null"
            );
        }
        const auto& endpoint = *static_cast<const NativeSdkEndpoint*>(context);
        if (!endpoint.phase_callbacks_admitted) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1,
                "model is not admitted for this phase callback"
            );
        }
        if ((std::isfinite(endpoint.temperature_min_k)
             && temperature_k < endpoint.temperature_min_k)
            || (std::isfinite(endpoint.temperature_max_k)
                && temperature_k > endpoint.temperature_max_k)) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1,
                "temperature is outside the parameter source domain"
            );
        }

        try {
            const auto block = endpoint.model.evaluate_pure_phase_parameters(
                temperature_k,
                amount_mol,
                volume_m3,
                segment_count,
                segment_diameter_angstrom,
                dispersion_energy_over_k
            );
            result->helmholtz_over_rt_reference_amount =
                block.helmholtz_over_rt_reference_amount;
            std::copy(block.gradient.begin(), block.gradient.end(), result->gradient);
            std::copy(block.hessian.begin(), block.hessian.end(), result->hessian);
            result->pressure_pa = block.pressure_pa;
            result->chemical_potential_over_rt = block.chemical_potential_over_rt;
            std::snprintf(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                "%s",
                endpoint.parameter_fingerprint.c_str()
            );
            result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
            return result->status;
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::logic_error& error) {
            return fail(EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1, error.what());
        } catch (const std::exception& error) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1, error.what());
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1,
                "unknown provider evaluation failure"
            );
        }
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_associating_parameters(
        const void* context,
        double temperature_k,
        double amount_mol,
        double volume_m3,
        double segment_count,
        double segment_diameter_angstrom,
        double dispersion_energy_over_k,
        double association_energy_over_k,
        double association_volume,
        epcsaft_associating_parameterized_phase_block_result_v1* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if (result->struct_size
            != sizeof(epcsaft_associating_parameterized_phase_block_result_v1)) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        const uint32_t result_size = result->struct_size;
        std::memset(result, 0, sizeof(*result));
        result->struct_size = result_size;
        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "model context is null"
            );
        }
        const std::array<double, 8> inputs{
            temperature_k,
            amount_mol,
            volume_m3,
            segment_count,
            segment_diameter_angstrom,
            dispersion_energy_over_k,
            association_energy_over_k,
            association_volume,
        };
        if (!std::all_of(inputs.begin(), inputs.end(), [](double value) {
                return std::isfinite(value) && value > 0.0;
            })) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "joint associating phase inputs must be positive and finite"
            );
        }
        const auto& endpoint = *static_cast<const NativeSdkEndpoint*>(context);
        if (!endpoint.phase_callbacks_admitted) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1,
                "model is not admitted for this phase callback"
            );
        }
        if ((std::isfinite(endpoint.temperature_min_k)
             && temperature_k < endpoint.temperature_min_k)
            || (std::isfinite(endpoint.temperature_max_k)
                && temperature_k > endpoint.temperature_max_k)) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1,
                "temperature is outside the parameter source domain"
            );
        }
        try {
            const auto block =
                endpoint.model
                    .evaluate_diameter_basis_associating_pure_phase_parameters(
                    temperature_k,
                    amount_mol,
                    volume_m3,
                    segment_count,
                    segment_diameter_angstrom,
                    dispersion_energy_over_k,
                    association_energy_over_k,
                    association_volume
                );
            result->helmholtz_over_rt_reference_amount =
                block.helmholtz_over_rt_reference_amount;
            std::copy(
                block.gradient.begin(), block.gradient.end(), result->gradient
            );
            std::copy(
                block.hessian.begin(), block.hessian.end(), result->hessian
            );
            result->pressure_pa = block.pressure_pa;
            result->chemical_potential_over_rt =
                block.chemical_potential_over_rt;
            std::snprintf(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                "%s",
                endpoint.parameter_fingerprint.c_str()
            );
            std::snprintf(
                result->topology_fingerprint,
                sizeof(result->topology_fingerprint),
                "%s",
                endpoint.topology_fingerprint.c_str()
            );
            result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
            return result->status;
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::logic_error& error) {
            return fail(EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1, error.what());
        } catch (const std::exception& error) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1, error.what());
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1,
                "unknown joint associating evaluation failure"
            );
        }
    }

    template <typename Evaluator>
    // EqExclude: transport and buffer mechanics
    static int evaluate_mixture_common(
        const void* context,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double volume_m3,
        std::size_t additional_coordinate_count,
        bool additional_inputs_valid,
        const char* invalid_state_message,
        epcsaft_mixture_phase_block_result_v1* result,
        Evaluator evaluator
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if (result->struct_size != sizeof(epcsaft_mixture_phase_block_result_v1)) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        if (context == nullptr) {
            result->status = EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
            std::snprintf(
                result->error,
                sizeof(result->error),
                "%s",
                "model context is null"
            );
            return result->status;
        }

        const auto& endpoint = *static_cast<const NativeSdkEndpoint*>(context);
        if (!endpoint.phase_callbacks_admitted) {
            result->status = EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1;
            std::snprintf(
                result->error,
                sizeof(result->error),
                "%s",
                "model is not admitted for this phase callback"
            );
            return result->status;
        }
        const std::size_t coordinate_count =
            endpoint.table.component_count + additional_coordinate_count;
        if (result->coordinate_count != coordinate_count
            || result->gradient_capacity < coordinate_count
            || result->hessian_capacity < coordinate_count * coordinate_count
            || result->gradient == nullptr || result->hessian == nullptr) {
            result->status = EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
            std::snprintf(
                result->error,
                sizeof(result->error),
                "%s",
                "mixture result buffers do not match the model coordinate count"
            );
            return result->status;
        }

        result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
        result->helmholtz_over_rt_reference_amount = 0.0;
        result->pressure_pa = 0.0;
        result->parameter_fingerprint[0] = '\0';
        result->error[0] = '\0';
        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };

        if (amounts_mol == nullptr || amount_count != endpoint.table.component_count) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "amount buffer does not match the model component count"
            );
        }
        if (!std::isfinite(temperature_k) || !std::isfinite(volume_m3)
            || temperature_k <= 0.0 || volume_m3 <= 0.0
            || !additional_inputs_valid
            || !std::all_of(amounts_mol, amounts_mol + amount_count, [](double amount) {
                return std::isfinite(amount) && amount > 0.0;
            })) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                invalid_state_message
            );
        }
        if ((std::isfinite(endpoint.temperature_min_k)
             && temperature_k < endpoint.temperature_min_k)
            || (std::isfinite(endpoint.temperature_max_k)
                && temperature_k > endpoint.temperature_max_k)) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1,
                "temperature is outside the parameter source domain"
            );
        }

        try {
            const std::vector<double> amounts(amounts_mol, amounts_mol + amount_count);
            const auto block = evaluator(endpoint.model, amounts);
            std::copy(block.gradient.begin(), block.gradient.end(), result->gradient);
            std::copy(block.hessian.begin(), block.hessian.end(), result->hessian);
            result->helmholtz_over_rt_reference_amount =
                block.helmholtz_over_rt_reference_amount;
            result->pressure_pa = block.pressure_pa;
            std::snprintf(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                "%s",
                endpoint.parameter_fingerprint.c_str()
            );
            result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
            return result->status;
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::logic_error& error) {
            return fail(EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1, error.what());
        } catch (const std::exception& error) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1, error.what());
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1,
                "unknown provider evaluation failure"
            );
        }
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_pure_parameter(
        const void* context,
        double temperature_k,
        double amount_mol,
        double volume_m3,
        uint32_t parameter_family,
        double parameter_value,
        epcsaft_mixture_phase_block_result_v1* result
    ) noexcept {
        std::optional<ActivePureParameter> family;
        switch (parameter_family) {
            case EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_COUNT_V1:
                family = ActivePureParameter::SegmentCount;
                break;
            case EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_DIAMETER_V1:
                family = ActivePureParameter::SegmentDiameter;
                break;
            case EPCSAFT_NATIVE_PARAMETER_FAMILY_DISPERSION_ENERGY_OVER_K_V1:
                family = ActivePureParameter::DispersionEnergyOverK;
                break;
            case EPCSAFT_NATIVE_PARAMETER_FAMILY_ASSOCIATION_ENERGY_OVER_K_V1:
                family = ActivePureParameter::AssociationEnergyOverK;
                break;
            case EPCSAFT_NATIVE_PARAMETER_FAMILY_ASSOCIATION_VOLUME_V1:
                family = ActivePureParameter::AssociationVolume;
                break;
            default:
                break;
        }
        const std::array<double, 1> amounts{amount_mol};
        return evaluate_mixture_common(
            context,
            temperature_k,
            amounts.data(),
            amounts.size(),
            volume_m3,
            2,
            family.has_value() && std::isfinite(parameter_value),
            "temperature, amount, volume, and active pure parameter must be "
            "finite and state coordinates positive",
            result,
            [temperature_k, amount_mol, volume_m3, family, parameter_value](
                const NativeModel& model,
                const std::vector<double>&
            ) {
                return model.evaluate_pure_phase_parameter(
                    temperature_k,
                    amount_mol,
                    volume_m3,
                    *family,
                    parameter_value
                );
            }
        );
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_mixture(
        const void* context,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double volume_m3,
        epcsaft_mixture_phase_block_result_v1* result
    ) noexcept {
        return evaluate_mixture_common(
            context,
            temperature_k,
            amounts_mol,
            amount_count,
            volume_m3,
            1,
            true,
            "temperature, amounts, and volume must be positive and finite",
            result,
            [temperature_k, volume_m3](
                const NativeModel& model,
                const std::vector<double>& amounts
            ) {
                return model.evaluate_mixture_phase(
                    temperature_k, amounts, volume_m3
                );
            }
        );
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_mixture_kij(
        const void* context,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double volume_m3,
        double k_ij,
        epcsaft_mixture_phase_block_result_v1* result
    ) noexcept {
        return evaluate_mixture_common(
            context,
            temperature_k,
            amounts_mol,
            amount_count,
            volume_m3,
            2,
            std::isfinite(k_ij),
            "temperature, amounts, volume, and k_ij must be finite and state coordinates positive",
            result,
            [temperature_k, volume_m3, k_ij](
                const NativeModel& model,
                const std::vector<double>& amounts
            ) {
                return model.evaluate_mixture_phase_kij(
                    temperature_k, amounts, volume_m3, k_ij
                );
            }
        );
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_mixture_lij(
        const void* context,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double volume_m3,
        double l_ij,
        epcsaft_mixture_phase_block_result_v1* result
    ) noexcept {
        return evaluate_mixture_common(
            context,
            temperature_k,
            amounts_mol,
            amount_count,
            volume_m3,
            2,
            std::isfinite(l_ij),
            "temperature, amounts, volume, and l_ij must be finite and state coordinates positive",
            result,
            [temperature_k, volume_m3, l_ij](
                const NativeModel& model,
                const std::vector<double>& amounts
            ) {
                return model.evaluate_mixture_phase_lij(
                    temperature_k, amounts, volume_m3, l_ij
                );
            }
        );
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_electrolyte(
        const void* context,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double volume_m3,
        epcsaft_mixture_phase_block_result_v1* result
    ) noexcept {
        const auto* endpoint = static_cast<const NativeSdkEndpoint*>(context);
        return evaluate_mixture_common(
            context,
            temperature_k,
            amounts_mol,
            amount_count,
            volume_m3,
            1,
            true,
            "temperature, amounts, and volume must be positive and finite",
            result,
            [temperature_k, volume_m3, endpoint](
                const NativeModel& model,
                const std::vector<double>& amounts
            ) {
                return model.evaluate_electrolyte_phase(
                    temperature_k,
                    amounts,
                    volume_m3,
                    endpoint->ion_mole_fraction_max
                );
            }
        );
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_molar_volume_bounds(
        const void* context,
        const char* expected_parameter_fingerprint,
        double temperature_k,
        const double* mole_fractions,
        size_t component_count,
        double packing_fraction_min,
        double packing_fraction_max,
        double* molar_volume_bounds_m3_per_mol,
        size_t bound_count
    ) noexcept {
        if (bound_count != EPCSAFT_NATIVE_MOLAR_VOLUME_BOUND_COUNT_V1) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        if (context == nullptr || expected_parameter_fingerprint == nullptr
            || mole_fractions == nullptr
            || molar_volume_bounds_m3_per_mol == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        const auto& endpoint = *static_cast<const NativeSdkEndpoint*>(context);
        if (!endpoint.phase_callbacks_admitted) {
            return EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1;
        }
        if (std::strcmp(
                expected_parameter_fingerprint,
                endpoint.parameter_fingerprint.c_str()
            ) != 0) {
            return EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1;
        }
        if (component_count != endpoint.table.component_count
            || !std::isfinite(temperature_k) || temperature_k <= 0.0
            || !std::isfinite(packing_fraction_min)
            || !std::isfinite(packing_fraction_max)
            || packing_fraction_min <= 0.0
            || packing_fraction_min >= packing_fraction_max) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if ((std::isfinite(endpoint.temperature_min_k)
             && temperature_k < endpoint.temperature_min_k)
            || (std::isfinite(endpoint.temperature_max_k)
                && temperature_k > endpoint.temperature_max_k)) {
            return EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1;
        }

        try {
            const std::vector<double> composition(
                mole_fractions, mole_fractions + component_count
            );
            const auto bounds = endpoint.model.molar_volume_bounds(
                temperature_k,
                composition,
                packing_fraction_min,
                packing_fraction_max,
                endpoint.ion_mole_fraction_max
            );
            std::copy(
                bounds.begin(), bounds.end(), molar_volume_bounds_m3_per_mol
            );
            return EPCSAFT_NATIVE_STATUS_OK_V1;
        } catch (const std::invalid_argument&) {
            return EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1;
        } catch (...) {
            return EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1;
        }
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_packing_fraction(
        const void* context,
        const char* expected_parameter_fingerprint,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double volume_m3,
        double* packing_fraction,
        double* gradient,
        size_t gradient_capacity,
        double* hessian,
        size_t hessian_capacity
    ) noexcept {
        if (context == nullptr || expected_parameter_fingerprint == nullptr
            || amounts_mol == nullptr || packing_fraction == nullptr
            || gradient == nullptr || hessian == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        const auto& endpoint = *static_cast<const NativeSdkEndpoint*>(context);
        if (!endpoint.phase_callbacks_admitted) {
            return EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1;
        }
        const std::size_t coordinate_count = endpoint.table.component_count + 1;
        if (gradient_capacity != coordinate_count
            || hessian_capacity != coordinate_count * coordinate_count) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        if (std::strcmp(
                expected_parameter_fingerprint,
                endpoint.parameter_fingerprint.c_str()
            ) != 0) {
            return EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1;
        }
        if (amount_count != endpoint.table.component_count
            || !std::isfinite(temperature_k) || temperature_k <= 0.0
            || !std::isfinite(volume_m3) || volume_m3 <= 0.0
            || !std::all_of(amounts_mol, amounts_mol + amount_count, [](double amount) {
                return std::isfinite(amount) && amount > 0.0;
            })) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if ((std::isfinite(endpoint.temperature_min_k)
             && temperature_k < endpoint.temperature_min_k)
            || (std::isfinite(endpoint.temperature_max_k)
                && temperature_k > endpoint.temperature_max_k)) {
            return EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1;
        }

        try {
            const std::vector<double> amounts(amounts_mol, amounts_mol + amount_count);
            const auto block = endpoint.model.evaluate_packing_fraction(
                temperature_k,
                amounts,
                volume_m3,
                endpoint.ion_mole_fraction_max
            );
            *packing_fraction = block.packing_fraction;
            std::copy(block.gradient.begin(), block.gradient.end(), gradient);
            std::copy(block.hessian.begin(), block.hessian.end(), hessian);
            return EPCSAFT_NATIVE_STATUS_OK_V1;
        } catch (const std::invalid_argument&) {
            return EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1;
        } catch (const std::logic_error&) {
            return EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1;
        } catch (...) {
            return EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1;
        }
    }
};

}  // namespace

// EqExclude: transport and buffer mechanics
const epcsaft_native_sdk_v1* create_native_sdk_v1(
    const NativeModel& model,
    const std::string& parameter_fingerprint,
    double temperature_min_k,
    double temperature_max_k,
    const std::vector<std::string>& component_ids,
    double ion_mole_fraction_max,
    bool phase_callbacks_admitted,
    bool active_born_admitted,
    bool aqueous_kij_admitted,
    bool ion_fraction_suppression_admitted,
    double pressure_min_pa,
    double pressure_max_pa,
    bool neutral_reference_admitted,
    bool unknown_applicability_admitted,
    bool aqueous_solvation_factor_admitted,
    bool ion_solvation_kij_admitted,
    bool staged_trial_bundle_admitted,
    const std::string& topology_fingerprint
) {
    if (!active_born_admitted && !aqueous_kij_admitted
        && !ion_fraction_suppression_admitted
        && !ion_solvation_kij_admitted
        && !neutral_reference_admitted && !aqueous_solvation_factor_admitted
        && (!model.supports_native_sdk_v1()
        && !model.supports_associating_pure_parameter_native_sdk_v1()
        && !model.supports_mixture_native_sdk_v1()
        && !model.supports_electrolyte_native_sdk_v1())) {
        throw std::invalid_argument(
            "native SDK v1 requires a supported phase-block model"
        );
    }
    if (active_born_admitted
        && !active_born_model_admitted(
            parameter_fingerprint,
            component_ids,
            model.charge_numbers(),
            staged_trial_bundle_admitted
        )) {
        throw std::invalid_argument("unsupported active Born model identity");
    }
    if (aqueous_kij_admitted
        && !aqueous_kij_model_admitted(
            parameter_fingerprint,
            component_ids,
            model.charge_numbers(),
            staged_trial_bundle_admitted
        )) {
        throw std::invalid_argument("unsupported aqueous active k_ij model identity");
    }
    if (ion_fraction_suppression_admitted
        && !ion_fraction_suppression_model_admitted(
            component_ids, model.charge_numbers()
        )) {
        throw std::invalid_argument(
            "unsupported active ion-fraction suppression model identity"
        );
    }
    if (ion_solvation_kij_admitted
        && !ion_solvation_kij_model_admitted(
            parameter_fingerprint,
            topology_fingerprint,
            component_ids, model.charge_numbers()
        )) {
        throw std::invalid_argument(
            "unsupported active ion-solvation k_ij model identity"
        );
    }
    if (aqueous_solvation_factor_admitted
        && !aqueous_solvation_factor_model_admitted(
            parameter_fingerprint,
            component_ids,
            model.charge_numbers(),
            staged_trial_bundle_admitted
        )) {
        throw std::invalid_argument(
            "unsupported aqueous active solvation-factor model identity"
        );
    }
    if (neutral_reference_admitted
        && !neutral_reference_model_shape_admitted(
            component_ids, model.charge_numbers()
        )) {
        throw std::invalid_argument("unsupported neutral reference model identity");
    }
    if (parameter_fingerprint.empty()
        || parameter_fingerprint.size() >= EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE) {
        throw std::invalid_argument("invalid parameter fingerprint");
    }
    if (!topology_fingerprint.empty()
        && (topology_fingerprint.size()
                >= EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE
            || topology_fingerprint.rfind("sha256:", 0) != 0
            || topology_fingerprint.size() != 71)) {
        throw std::invalid_argument("invalid topology fingerprint");
    }
    if (component_ids.size() != model.component_count()
        || std::any_of(
            component_ids.begin(), component_ids.end(), [](const auto& identifier) {
                return identifier.empty();
            }
        )) {
        throw std::invalid_argument("invalid component identifiers");
    }
    for (std::size_t left = 0; left < component_ids.size(); ++left) {
        for (std::size_t right = left + 1; right < component_ids.size(); ++right) {
            if (component_ids[left] == component_ids[right]) {
                throw std::invalid_argument("duplicate component identifier");
            }
        }
    }
    if (model.supports_electrolyte_native_sdk_v1()
        && ((!unknown_applicability_admitted
             && (std::isnan(temperature_min_k) || std::isnan(temperature_max_k)))
            || (!std::isnan(temperature_min_k)
             && (!std::isfinite(temperature_min_k) || temperature_min_k <= 0.0))
            || (!std::isnan(temperature_max_k)
                && (!std::isfinite(temperature_max_k) || temperature_max_k <= 0.0))
            || (std::isfinite(temperature_min_k)
                && std::isfinite(temperature_max_k)
                && temperature_min_k > temperature_max_k)
            || (!std::isnan(ion_mole_fraction_max)
                && (!std::isfinite(ion_mole_fraction_max)
                    || ion_mole_fraction_max <= 0.0
                    || ion_mole_fraction_max > 1.0)))) {
        throw std::invalid_argument("invalid electrolyte source applicability");
    }
    if (neutral_reference_admitted
        && ((!unknown_applicability_admitted
             && (std::isnan(pressure_min_pa) || std::isnan(pressure_max_pa)))
            || (!std::isnan(pressure_min_pa)
                && (!std::isfinite(pressure_min_pa) || pressure_min_pa <= 0.0))
            || (!std::isnan(pressure_max_pa)
                && (!std::isfinite(pressure_max_pa) || pressure_max_pa <= 0.0))
            || (std::isfinite(pressure_min_pa) && std::isfinite(pressure_max_pa)
            && pressure_min_pa > pressure_max_pa))) {
        throw std::invalid_argument("invalid neutral reference pressure applicability");
    }
    auto* endpoint = new NativeSdkEndpoint(
        model,
        parameter_fingerprint,
        temperature_min_k,
        temperature_max_k,
        component_ids,
        ion_mole_fraction_max,
        phase_callbacks_admitted,
        active_born_admitted,
        aqueous_kij_admitted,
        ion_fraction_suppression_admitted,
        ion_solvation_kij_admitted,
        pressure_min_pa,
        pressure_max_pa,
        neutral_reference_admitted,
        aqueous_solvation_factor_admitted,
        topology_fingerprint
    );
    return &endpoint->table;
}

// EqExclude: transport and buffer mechanics
void destroy_native_sdk_v1(const epcsaft_native_sdk_v1* sdk) noexcept {
    if (sdk != nullptr) {
        delete static_cast<const NativeSdkEndpoint*>(sdk->model_context);
    }
}

}  // namespace epcsaft
