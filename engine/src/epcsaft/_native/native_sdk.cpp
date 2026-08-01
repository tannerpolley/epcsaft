#include "native_sdk.hpp"
#include "density.hpp"
#include "electrolyte_reference.hpp"
#include "phase_block.hpp"

#include <epcsaft/native_sdk_v1.h>

#include <algorithm>
#include <array>
#include <cctype>
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
    const std::vector<std::string>& component_ids,
    const std::vector<int>& charges
) {
    return component_ids.size() == 3 && charges.size() == 3
        && charges[0] == 0 && charges[1] != 0 && charges[2] == -charges[1];
}

// EqExclude: input and domain validation
bool aqueous_kij_model_admitted(
    const std::vector<std::string>& component_ids,
    const std::vector<int>& charges
) {
    return component_ids.size() == 3
        && charges == std::vector<int>({0, 1, -1});
}

// EqExclude: input and domain validation
bool ion_fraction_suppression_model_admitted(
    const std::vector<std::string>& component_ids,
    const std::vector<int>& charges
) {
    return component_ids.size() == 3
        && charges == std::vector<int>({0, 1, -1});
}

// EqExclude: input and domain validation
bool ion_solvation_kij_model_admitted(
    const std::vector<std::string>& component_ids,
    const std::vector<int>& charges
) {
    return component_ids.size() == 3
        && charges == std::vector<int>({0, 1, -1});
}

// EqExclude: input and domain validation
bool aqueous_solvation_factor_model_admitted(
    const std::vector<std::string>& component_ids,
    const std::vector<int>& charges
) {
    return component_ids.size() == 3
        && charges == std::vector<int>({0, 1, -1});
}

// EqExclude: input and domain validation
struct NeutralReferenceBasis final {
    std::vector<std::size_t> neutral_components;
    std::vector<std::size_t> cations;
    std::vector<std::size_t> anions;
    std::vector<double> reference_composition;
    std::vector<double> rows;

    [[nodiscard]] std::size_t row_count() const {
        return reference_composition.empty()
            ? 0u
            : rows.size() / reference_composition.size();
    }
};

// EqExclude: input and domain validation
bool full_row_rank(
    const std::vector<double>& values,
    std::size_t rows,
    std::size_t columns
) {
    if (rows == 0 || columns == 0 || values.size() != rows * columns) {
        return false;
    }
    std::vector<double> matrix(values);
    std::size_t rank = 0;
    for (std::size_t column = 0; column < columns && rank < rows; ++column) {
        std::size_t pivot = rank;
        while (pivot < rows
               && std::abs(matrix[pivot * columns + column]) <= 1.0e-12) {
            ++pivot;
        }
        if (pivot == rows) {
            continue;
        }
        for (std::size_t item = column; item < columns; ++item) {
            std::swap(
                matrix[rank * columns + item],
                matrix[pivot * columns + item]
            );
        }
        const double divisor = matrix[rank * columns + column];
        for (std::size_t item = column; item < columns; ++item) {
            matrix[rank * columns + item] /= divisor;
        }
        for (std::size_t row = 0; row < rows; ++row) {
            if (row == rank) {
                continue;
            }
            const double factor = matrix[row * columns + column];
            for (std::size_t item = column; item < columns; ++item) {
                matrix[row * columns + item] -=
                    factor * matrix[rank * columns + item];
            }
        }
        ++rank;
    }
    return rank == rows;
}

// EqID: neutral_reference_stoichiometric_basis
NeutralReferenceBasis make_neutral_reference_basis(
    const std::vector<std::string>& component_ids,
    const std::vector<int>& charges
) {
    NeutralReferenceBasis result;
    const std::size_t components = component_ids.size();
    if (components < 3 || components != charges.size()
        || std::any_of(
            component_ids.begin(), component_ids.end(), [](const auto& value) {
                return value.empty();
            }
        )
        || std::any_of(charges.begin(), charges.end(), [](int charge) {
            return charge == std::numeric_limits<int>::min();
        })) {
        return result;
    }
    std::vector<std::size_t> order(components);
    std::iota(order.begin(), order.end(), 0u);
    std::sort(order.begin(), order.end(), [&](std::size_t left, std::size_t right) {
        return component_ids[left] < component_ids[right];
    });
    for (std::size_t index = 1; index < order.size(); ++index) {
        if (component_ids[order[index - 1]] == component_ids[order[index]]) {
            return {};
        }
    }
    for (std::size_t component : order) {
        if (charges[component] == 0) {
            result.neutral_components.push_back(component);
        } else if (charges[component] > 0) {
            result.cations.push_back(component);
        } else {
            result.anions.push_back(component);
        }
    }
    if (result.neutral_components.empty() || result.cations.empty()
        || result.anions.empty()) {
        return {};
    }
    result.reference_composition.assign(components, 0.0);
    const double neutral_fraction =
        1.0 / static_cast<double>(result.neutral_components.size());
    for (std::size_t component : result.neutral_components) {
        result.reference_composition[component] = neutral_fraction;
    }
    auto append_row = [&](std::size_t left, double left_value,
                          std::size_t right = 0u, double right_value = 0.0) {
        result.rows.resize(result.rows.size() + components, 0.0);
        const std::size_t offset = result.rows.size() - components;
        result.rows[offset + left] = left_value;
        if (right_value != 0.0) {
            result.rows[offset + right] = right_value;
        }
    };
    for (std::size_t component : result.neutral_components) {
        append_row(component, 1.0);
    }
    auto append_pair = [&](std::size_t cation, std::size_t anion) {
        const int divisor = std::gcd(
            std::abs(charges[cation]), std::abs(charges[anion])
        );
        append_row(
            cation,
            static_cast<double>(std::abs(charges[anion]) / divisor),
            anion,
            static_cast<double>(std::abs(charges[cation]) / divisor)
        );
    };
    for (std::size_t anion : result.anions) {
        append_pair(result.cations.front(), anion);
    }
    for (std::size_t index = 1; index < result.cations.size(); ++index) {
        append_pair(result.cations[index], result.anions.front());
    }
    const std::size_t expected_rows = components - 1u;
    if (result.row_count() != expected_rows
        || !full_row_rank(result.rows, expected_rows, components)) {
        return {};
    }
    for (std::size_t row = 0; row < expected_rows; ++row) {
        double net_charge = 0.0;
        for (std::size_t component = 0; component < components; ++component) {
            net_charge += result.rows[row * components + component]
                * static_cast<double>(charges[component]);
        }
        if (net_charge != 0.0) {
            return {};
        }
    }
    return result;
}

// EqExclude: input and domain validation
bool neutral_reference_model_shape_admitted(
    const std::vector<std::string>& component_ids,
    const std::vector<int>& charges
) {
    return !component_ids.empty()
        && make_neutral_reference_basis(component_ids, charges).row_count()
            == component_ids.size() - 1u;
}

std::vector<std::pair<std::size_t, std::size_t>>
fixed_topology_association_pairs(const NativeModel& model) {
    const std::size_t sites = model.site_components().size();
    const auto& active = model.association_active_pairs();
    std::vector<std::pair<std::size_t, std::size_t>> pairs;
    for (std::size_t left = 0; left < sites; ++left) {
        for (std::size_t right = left; right < sites; ++right) {
            const std::size_t forward = left * sites + right;
            const std::size_t reverse = right * sites + left;
            if (active[forward] != 0 && active[reverse] != 0) {
                pairs.emplace_back(left, right);
            }
        }
    }
    return pairs;
}

// EqExclude: input and domain validation
bool fixed_topology_site_ids_valid(
    const NativeModel& model,
    const std::vector<std::string>& site_ids
) {
    if (!model.supports_fixed_topology_association_native_sdk_v1()
        || site_ids.size() != model.site_components().size()
        || std::any_of(
            site_ids.begin(), site_ids.end(), [](const auto& identifier) {
                return identifier.empty();
            }
        )) {
        return false;
    }
    for (std::size_t left = 0; left < site_ids.size(); ++left) {
        for (std::size_t right = left + 1; right < site_ids.size(); ++right) {
            if (model.site_components()[left] == model.site_components()[right]
                && site_ids[left] == site_ids[right]) {
                return false;
            }
        }
    }
    return true;
}

// EqExclude: input and domain validation
bool sha256_fingerprint_valid(const std::string& fingerprint) {
    return fingerprint.size() == 71
        && fingerprint.rfind("sha256:", 0) == 0
        && std::all_of(
            fingerprint.begin() + 7,
            fingerprint.end(),
            [](unsigned char character) {
                return std::isxdigit(character) != 0;
            }
        );
}

// EqExclude: transport and buffer mechanics
std::vector<std::string> make_fixed_topology_parameter_ids(
    const NativeModel& model,
    const std::vector<std::string>& component_ids,
    const std::vector<std::string>& site_ids
) {
    if (!fixed_topology_site_ids_valid(model, site_ids)) {
        return {};
    }
    std::vector<std::string> result{
        "component:" + component_ids[0] + "/segment_count",
        "component:" + component_ids[0] + "/segment_diameter",
        "component:" + component_ids[0] + "/dispersion_energy_over_k",
    };
    const auto& site_components = model.site_components();
    for (const auto& [left, right] : fixed_topology_association_pairs(model)) {
        const std::string endpoints =
            "association:"
            + component_ids[static_cast<std::size_t>(site_components[left])]
            + "/" + site_ids[left] + "--"
            + component_ids[static_cast<std::size_t>(site_components[right])]
            + "/" + site_ids[right] + "/";
        result.push_back(endpoints + "association_energy_over_k");
        result.push_back(endpoints + "association_volume");
    }
    return result;
}

std::vector<epcsaft_fixed_topology_parameter_slot_v1>
make_fixed_topology_parameter_slots(
    const NativeModel& model,
    const std::vector<std::string>& parameter_ids,
    const std::vector<std::string>& site_ids
) {
    if (!fixed_topology_site_ids_valid(model, site_ids)) {
        return {};
    }
    std::vector<epcsaft_fixed_topology_parameter_slot_v1> result;
    result.reserve(parameter_ids.size());
    const auto add_component = [&](uint32_t family, const char* unit) {
        const std::size_t index = result.size();
        result.push_back({
            sizeof(epcsaft_fixed_topology_parameter_slot_v1),
            family,
            0,
            -1,
            -1,
            -1,
            -1,
            parameter_ids[index].c_str(),
            unit,
            0.0,
            std::numeric_limits<double>::infinity(),
        });
    };
    add_component(
        EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_COUNT_V1,
        "dimensionless"
    );
    add_component(
        EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_DIAMETER_V1,
        "angstrom"
    );
    add_component(
        EPCSAFT_NATIVE_PARAMETER_FAMILY_DISPERSION_ENERGY_OVER_K_V1,
        "kelvin"
    );
    const auto& site_components = model.site_components();
    for (const auto& [left, right] : fixed_topology_association_pairs(model)) {
        const auto add_association = [&](uint32_t family, const char* unit) {
            const std::size_t index = result.size();
            result.push_back({
                sizeof(epcsaft_fixed_topology_parameter_slot_v1),
                family,
                -1,
                static_cast<int32_t>(left),
                static_cast<int32_t>(right),
                static_cast<int32_t>(site_components[left]),
                static_cast<int32_t>(site_components[right]),
                parameter_ids[index].c_str(),
                unit,
                0.0,
                std::numeric_limits<double>::infinity(),
            });
        };
        add_association(
            EPCSAFT_NATIVE_PARAMETER_FAMILY_ASSOCIATION_ENERGY_OVER_K_V1,
            "kelvin"
        );
        add_association(
            EPCSAFT_NATIVE_PARAMETER_FAMILY_ASSOCIATION_VOLUME_V1,
            "dimensionless"
        );
    }
    return result;
}

struct NativeSdkEndpoint final {
    NativeModel model;
    std::string parameter_fingerprint;
    std::string topology_fingerprint;
    std::string artifact_fingerprint;
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
    NeutralReferenceBasis neutral_reference_basis;
    std::vector<std::string> site_ids;
    std::vector<const char*> site_id_pointers;
    std::vector<int32_t> site_components;
    std::vector<int32_t> site_multiplicities;
    std::vector<std::string> fixed_topology_parameter_ids;
    std::vector<epcsaft_fixed_topology_parameter_slot_v1>
        fixed_topology_parameter_slots;
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
    std::vector<std::array<epcsaft_native_capability_coordinate_v1, 1>>
        reacting_parameter_coordinates;
    std::vector<epcsaft_native_capability_descriptor_v1> capabilities;
    epcsaft_fixed_topology_association_descriptor_v1
        fixed_topology_association;
    epcsaft_native_sdk_v1 table;

    NativeSdkEndpoint(
        const NativeModel& source_model,
        std::string fingerprint,
        double minimum_temperature,
        double maximum_temperature,
        std::vector<std::string> identifiers,
        std::vector<std::string> site_identifiers,
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
        std::string topology_identity,
        std::string artifact_identity
    )
        : model(source_model),
          parameter_fingerprint(std::move(fingerprint)),
          topology_fingerprint(std::move(topology_identity)),
          artifact_fingerprint(std::move(artifact_identity)),
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
          neutral_reference_basis(make_neutral_reference_basis(
              component_ids, source_model.charge_numbers()
          )),
          site_ids(std::move(site_identifiers)),
          site_id_pointers([&] {
              std::vector<const char*> pointers;
              pointers.reserve(site_ids.size());
              for (const auto& identifier : site_ids) {
                  pointers.push_back(identifier.c_str());
              }
              return pointers;
          }()),
          site_components(
              source_model.site_components().begin(),
              source_model.site_components().end()
          ),
          site_multiplicities(
              source_model.site_multiplicities().begin(),
              source_model.site_multiplicities().end()
          ),
          fixed_topology_parameter_ids(make_fixed_topology_parameter_ids(
              source_model, component_ids, site_ids
          )),
          fixed_topology_parameter_slots(make_fixed_topology_parameter_slots(
              source_model, fixed_topology_parameter_ids, site_ids
          )),
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
          reacting_parameter_coordinates([&] {
              std::vector<
                  std::array<epcsaft_native_capability_coordinate_v1, 1>
              > coordinates;
              coordinates.reserve(component_ids.size() * 2);
              for (std::size_t component = 0;
                   component < component_ids.size();
                   ++component) {
                  coordinates.push_back({{
                      {
                          sizeof(epcsaft_native_capability_coordinate_v1),
                          EPCSAFT_NATIVE_CAPABILITY_COORDINATE_SEGMENT_DIAMETER_V1,
                          static_cast<int32_t>(component),
                          -1,
                          -1,
                          "angstrom",
                      },
                  }});
                  coordinates.push_back({{
                      {
                          sizeof(epcsaft_native_capability_coordinate_v1),
                          EPCSAFT_NATIVE_CAPABILITY_COORDINATE_DISPERSION_ENERGY_OVER_K_V1,
                          static_cast<int32_t>(component),
                          -1,
                          -1,
                          "kelvin",
                      },
                  }});
              }
              return coordinates;
          }()),
          capabilities([&] {
              std::vector<epcsaft_native_capability_descriptor_v1> result;
              pure_segment_diameter_coordinates[2].kind =
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_SEGMENT_DIAMETER_V1;
              pure_segment_diameter_coordinates[2].unit = "angstrom";
              pure_dispersion_energy_coordinates[2].kind =
                  EPCSAFT_NATIVE_CAPABILITY_COORDINATE_DISPERSION_ENERGY_OVER_K_V1;
              pure_dispersion_energy_coordinates[2].unit = "kelvin";
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
              if (phase_callbacks_admitted
                  && source_model.supports_electrolyte_native_sdk_v1()
                  && !topology_fingerprint.empty()) {
                  for (std::size_t component = 0;
                       component < component_ids.size();
                       ++component) {
                      if (
                          source_model
                              .has_temperature_independent_segment_diameter(
                                  component
                              )
                      ) {
                          const auto& diameter =
                              reacting_parameter_coordinates[2 * component];
                          result.push_back(descriptor(
                              EPCSAFT_NATIVE_CAPABILITY_REACTING_COMPONENT_DIAMETER_HELMHOLTZ_V1,
                              EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_DIAMETER_V1,
                              EPCSAFT_NATIVE_PARAMETER_IDENTITY_COMPONENT_V1,
                              EPCSAFT_NATIVE_MODEL_DOMAIN_REACTING_ELECTROLYTE_PHASE_V1,
                              component_ids.size() + 1,
                              diameter
                          ));
                      }
                      const auto& energy =
                          reacting_parameter_coordinates[2 * component + 1];
                      result.push_back(descriptor(
                          EPCSAFT_NATIVE_CAPABILITY_REACTING_COMPONENT_DISPERSION_HELMHOLTZ_V1,
                          EPCSAFT_NATIVE_PARAMETER_FAMILY_DISPERSION_ENERGY_OVER_K_V1,
                          EPCSAFT_NATIVE_PARAMETER_IDENTITY_COMPONENT_V1,
                          EPCSAFT_NATIVE_MODEL_DOMAIN_REACTING_ELECTROLYTE_PHASE_V1,
                          component_ids.size() + 1,
                          energy
                      ));
                  }
              }
              return result;
          }()),
          fixed_topology_association([&] {
              epcsaft_fixed_topology_association_descriptor_v1 value{};
              if (!phase_callbacks_admitted
                  || !source_model
                          .supports_fixed_topology_association_native_sdk_v1()
                  || !fixed_topology_site_ids_valid(source_model, site_ids)
                  || !sha256_fingerprint_valid(artifact_fingerprint)
                  || topology_fingerprint.empty()) {
                  return value;
              }
              value.struct_size =
                  sizeof(epcsaft_fixed_topology_association_descriptor_v1);
              value.schema_version =
                  EPCSAFT_FIXED_TOPOLOGY_ASSOCIATION_SCHEMA_VERSION_V1;
              value.observation_contract =
                  EPCSAFT_NATIVE_OBSERVATION_FIXED_COMPOSITION_HELMHOLTZ_PHASE_V1;
              value.model_domain =
                  EPCSAFT_NATIVE_MODEL_DOMAIN_NEUTRAL_ASSOCIATING_PURE_V1;
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
              value.state_coordinate_count = 2;
              value.slot_count = fixed_topology_parameter_slots.size();
              value.slots = fixed_topology_parameter_slots.data();
              value.component_count = component_ids.size();
              value.component_ids = component_id_pointers.data();
              value.site_count = site_ids.size();
              value.site_ids = site_id_pointers.data();
              value.site_component_indices = site_components.data();
              value.site_multiplicities = site_multiplicities.data();
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
                  value.artifact_fingerprint,
                  sizeof(value.artifact_fingerprint),
                  "%s",
                  artifact_fingerprint.c_str()
              );
              std::snprintf(
                  value.helmholtz_basis_id,
                  sizeof(value.helmholtz_basis_id),
                  "%s",
                  EPCSAFT_NATIVE_SDK_V1_HELMHOLTZ_BASIS_ID
              );
              std::snprintf(
                  value.association_basis_id,
                  sizeof(value.association_basis_id),
                  "%s",
                  "ordinary_sigma_cubed_kappa"
              );
              return value;
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
              neutral_reference_admitted
                  ? neutral_reference_basis.row_count()
                  : 0u,
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
              0u,
              nullptr,
              sizeof(epcsaft_electrolyte_phase_value_result_v1),
              &NativeSdkEndpoint::evaluate_electrolyte_value,
              sizeof(epcsaft_neutral_reference_derivative_result_v1),
              &NativeSdkEndpoint::evaluate_neutral_reference_derivatives,
              sizeof(epcsaft_reacting_phase_parameter_result_v1),
              &NativeSdkEndpoint::evaluate_reacting_phase_parameters,
              sizeof(epcsaft_inverse_packing_geometry_result_v1),
              &NativeSdkEndpoint::evaluate_inverse_packing_geometry,
              nullptr,
              sizeof(epcsaft_fixed_topology_phase_result_v1),
              fixed_topology_association.struct_size == 0
                  ? nullptr
                  : &fixed_topology_association,
              &NativeSdkEndpoint::evaluate_fixed_topology_association,
          } {}

    // EqExclude: transport and buffer mechanics
    static std::vector<NativeActiveParameter> active_parameters(
        const NativeSdkEndpoint& endpoint,
        const epcsaft_active_parameter_request_v1* requests,
        std::size_t request_count
    ) {
        if (request_count == 0) {
            return {};
        }
        if (requests == nullptr) {
            throw std::invalid_argument(
                "active parameter request pointer is null"
            );
        }
        std::vector<NativeActiveParameter> result;
        result.reserve(request_count);
        for (std::size_t index = 0; index < request_count; ++index) {
            const auto& request = requests[index];
            if (request.struct_size
                    != sizeof(epcsaft_active_parameter_request_v1)
                || request.parameter_identity
                    != EPCSAFT_NATIVE_PARAMETER_IDENTITY_COMPONENT_V1
                || request.component_index < 0
                || static_cast<std::size_t>(request.component_index)
                    >= endpoint.component_ids.size()
                || request.pair_component_index_a != -1
                || request.pair_component_index_b != -1
                || request.unit == nullptr || !std::isfinite(request.value)
                || request.value <= 0.0) {
                throw std::invalid_argument(
                    "active parameter request is malformed"
                );
            }
            ActiveComponentParameter family;
            const char* expected_unit = nullptr;
            uint32_t expected_capability = 0;
            uint32_t expected_coordinate_kind = 0;
            if (request.parameter_family
                == EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_DIAMETER_V1) {
                family = ActiveComponentParameter::SegmentDiameter;
                expected_unit = "angstrom";
                expected_capability =
                    EPCSAFT_NATIVE_CAPABILITY_REACTING_COMPONENT_DIAMETER_HELMHOLTZ_V1;
                expected_coordinate_kind =
                    EPCSAFT_NATIVE_CAPABILITY_COORDINATE_SEGMENT_DIAMETER_V1;
            } else if (
                request.parameter_family
                == EPCSAFT_NATIVE_PARAMETER_FAMILY_DISPERSION_ENERGY_OVER_K_V1
            ) {
                family = ActiveComponentParameter::DispersionEnergyOverK;
                expected_unit = "kelvin";
                expected_capability =
                    EPCSAFT_NATIVE_CAPABILITY_REACTING_COMPONENT_DISPERSION_HELMHOLTZ_V1;
                expected_coordinate_kind =
                    EPCSAFT_NATIVE_CAPABILITY_COORDINATE_DISPERSION_ENERGY_OVER_K_V1;
            } else {
                throw std::logic_error(
                    "active parameter family is not advertised"
                );
            }
            if (std::strcmp(request.unit, expected_unit) != 0) {
                throw std::invalid_argument(
                    "active parameter unit does not match its descriptor"
                );
            }
            const bool advertised = std::any_of(
                endpoint.capabilities.begin(),
                endpoint.capabilities.end(),
                [&](const epcsaft_native_capability_descriptor_v1& descriptor) {
                    if (descriptor.capability != expected_capability
                        || descriptor.parameter_family
                            != request.parameter_family
                        || descriptor.parameter_identity
                            != request.parameter_identity
                        || descriptor.model_domain
                            != EPCSAFT_NATIVE_MODEL_DOMAIN_REACTING_ELECTROLYTE_PHASE_V1
                        || descriptor.active_parameter_count != 1
                        || descriptor.coordinate_count != 1
                        || descriptor.coordinates == nullptr
                        || std::strcmp(
                            descriptor.parameter_fingerprint,
                            endpoint.parameter_fingerprint.c_str()
                        ) != 0
                        || std::strcmp(
                            descriptor.topology_fingerprint,
                            endpoint.topology_fingerprint.c_str()
                        ) != 0) {
                        return false;
                    }
                    const auto& coordinate = descriptor.coordinates[0];
                    return coordinate.kind == expected_coordinate_kind
                        && coordinate.component_index
                            == request.component_index
                        && coordinate.pair_component_index_a
                            == request.pair_component_index_a
                        && coordinate.pair_component_index_b
                            == request.pair_component_index_b
                        && coordinate.unit != nullptr
                        && std::strcmp(coordinate.unit, request.unit) == 0;
                }
            );
            if (!advertised) {
                throw std::logic_error(
                    "active parameter request is not advertised by this model"
                );
            }
            if (std::any_of(
                    result.begin(),
                    result.end(),
                    [&](const NativeActiveParameter& item) {
                        return item.family == family
                            && item.component
                                == static_cast<std::size_t>(
                                    request.component_index
                                );
                    }
                )) {
                throw std::invalid_argument(
                    "active parameter request is duplicated"
                );
            }
            result.push_back({
                family,
                static_cast<std::size_t>(request.component_index),
                request.value,
            });
        }
        return result;
    }

    // EqID: phase_packing_fraction
    // EqExclude: transport and buffer mechanics
    static int evaluate_inverse_packing_geometry(
        const void* context,
        const char* expected_fingerprint,
        const char* expected_topology_fingerprint,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double log_packing_fraction,
        const epcsaft_active_parameter_request_v1* requests,
        size_t request_count,
        epcsaft_inverse_packing_geometry_result_v1* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if (result->struct_size
            != sizeof(epcsaft_inverse_packing_geometry_result_v1)) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
        result->volume_m3 = std::numeric_limits<double>::quiet_NaN();
        result->parameter_fingerprint[0] = '\0';
        result->topology_fingerprint[0] = '\0';
        result->error[0] = '\0';
        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr || expected_fingerprint == nullptr
            || expected_topology_fingerprint == nullptr
            || amounts_mol == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "inverse packing required pointer is null"
            );
        }
        const auto& endpoint =
            *static_cast<const NativeSdkEndpoint*>(context);
        if (endpoint.topology_fingerprint.empty()
            || std::strcmp(
                expected_fingerprint,
                endpoint.parameter_fingerprint.c_str()
            ) != 0
            || std::strcmp(
                expected_topology_fingerprint,
                endpoint.topology_fingerprint.c_str()
            ) != 0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1,
                "inverse packing model or topology is not admitted"
            );
        }
        if (request_count > 2 * endpoint.table.component_count) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "inverse packing active parameter count exceeds model maximum"
            );
        }
        const std::size_t coordinate_count =
            endpoint.table.component_count + 1 + request_count;
        if (amount_count != endpoint.table.component_count
            || result->coordinate_count != coordinate_count
            || result->active_parameter_count != request_count
            || result->gradient_capacity != coordinate_count
            || result->hessian_capacity
                != coordinate_count * coordinate_count
            || result->gradient == nullptr || result->hessian == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1,
                "inverse packing shape or capacity does not match"
            );
        }
        if (!std::isfinite(temperature_k) || temperature_k <= 0.0
            || !std::all_of(
                amounts_mol,
                amounts_mol + amount_count,
                [](double amount) {
                    return std::isfinite(amount) && amount > 0.0;
                }
            )
            || !std::isfinite(log_packing_fraction)) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "inverse packing state is invalid"
            );
        }
        const double packing_fraction = std::exp(log_packing_fraction);
        if (!std::isfinite(packing_fraction) || packing_fraction <= 0.0
            || packing_fraction >= detail::close_packing_fraction()) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1,
                "inverse packing fraction is outside the model domain"
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
            const auto parameters =
                active_parameters(endpoint, requests, request_count);
            const std::vector<double> amounts(
                amounts_mol, amounts_mol + amount_count
            );
            const auto block =
                endpoint.model.evaluate_inverse_packing_geometry(
                    temperature_k,
                    amounts,
                    log_packing_fraction,
                    parameters
                );
            std::copy(
                block.gradient.begin(),
                block.gradient.end(),
                result->gradient
            );
            std::copy(
                block.hessian.begin(),
                block.hessian.end(),
                result->hessian
            );
            result->volume_m3 = block.volume_m3;
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
            return result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::logic_error& error) {
            return fail(EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1, error.what());
        } catch (const std::exception& error) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1, error.what());
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1,
                "unknown inverse packing evaluation failure"
            );
        }
    }

    // EqID: neutral_reference_stoichiometric_basis
    // EqID: fixed_pressure_total_derivative
    static int evaluate_neutral_reference_derivatives(
        const void* context,
        const char* expected_fingerprint,
        const char* expected_topology_fingerprint,
        double temperature_k,
        double pressure_pa,
        const epcsaft_active_parameter_request_v1* requests,
        size_t request_count,
        epcsaft_neutral_reference_derivative_result_v1* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if (result->struct_size
            != sizeof(epcsaft_neutral_reference_derivative_result_v1)) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        const double nan = std::numeric_limits<double>::quiet_NaN();
        result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
        result->derivative_availability =
            EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_NONE_V1;
        result->source_pressure_min_pa = nan;
        result->source_pressure_max_pa = nan;
        result->maximum_root_residual_pa = nan;
        result->minimum_pressure_density_derivative_pa_m3_per_mol = nan;
        result->maximum_density_condition_number = nan;
        result->reference_derivative_convergence_error = nan;
        result->maximum_relative_root_bracket_width = nan;
        result->maximum_relative_root_density_step = nan;
        result->stable_root_count = 0;
        result->selected_stable_root_index = 0;
        result->reference_branch[0] = '\0';
        result->parameter_fingerprint[0] = '\0';
        result->topology_fingerprint[0] = '\0';
        result->helmholtz_basis_id[0] = '\0';
        result->error[0] = '\0';
        auto fail = [&](int status, const char* message) {
            result->status = status;
            result->value.status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr || expected_fingerprint == nullptr
            || expected_topology_fingerprint == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "neutral-reference derivative required pointer is null"
            );
        }
        const auto& endpoint =
            *static_cast<const NativeSdkEndpoint*>(context);
        if (!endpoint.neutral_reference_admitted
            || endpoint.topology_fingerprint.empty()
            || std::strcmp(
                expected_fingerprint,
                endpoint.parameter_fingerprint.c_str()
            ) != 0
            || std::strcmp(
                expected_topology_fingerprint,
                endpoint.topology_fingerprint.c_str()
            ) != 0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1,
                "neutral-reference derivative model identity is unsupported"
            );
        }
        if (!std::isfinite(temperature_k) || temperature_k <= 0.0
            || !std::isfinite(pressure_pa) || pressure_pa <= 0.0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "temperature and pressure must be positive and finite"
            );
        }
        if (temperature_k < endpoint.temperature_min_k
            || temperature_k > endpoint.temperature_max_k
            || pressure_pa < endpoint.pressure_min_pa
            || pressure_pa > endpoint.pressure_max_pa) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1,
                "neutral-reference derivative state is outside the source domain"
            );
        }
        const std::size_t component_count = endpoint.component_ids.size();
        const auto& reference_basis = endpoint.neutral_reference_basis;
        const std::size_t row_count = reference_basis.row_count();
        auto& value = result->value;
        if (value.struct_size != sizeof(epcsaft_neutral_reference_result_v1)
            || value.component_count != component_count
            || value.neutral_basis_row_count != row_count
            || value.neutral_basis_capacity != row_count * component_count
            || value.contraction_capacity != row_count
            || value.reference_composition_capacity != component_count
            || value.neutral_basis == nullptr
            || value.log_fugacity_contractions == nullptr
            || value.reference_composition == nullptr
            || result->active_parameter_count != request_count
            || result->pressure_derivative_capacity != row_count
            || result->parameter_derivative_capacity
                != row_count * request_count
            || result->pressure_derivatives_per_pa == nullptr
            || (request_count != 0
                && result->parameter_derivatives == nullptr)) {
            return fail(
                EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1,
                "neutral-reference derivative buffers do not match immutable dimensions"
            );
        }
        try {
            const auto parameters = active_parameters(
                endpoint, requests, request_count
            );
            if (row_count != component_count - 1u) {
                throw std::logic_error(
                    "neutral-reference derivative topology is incomplete"
                );
            }
            const std::vector<double> reference_composition =
                reference_basis.reference_composition;
            std::vector<double> basis = reference_basis.rows;
            std::vector<double> contractions(row_count);
            std::vector<double> pressure_derivatives(row_count);
            std::vector<double> parameter_derivatives(
                row_count * request_count
            );
            double maximum_root_residual = 0.0;
            double minimum_slope = std::numeric_limits<double>::infinity();
            double maximum_condition = 0.0;
            double maximum_relative_bracket_width = 0.0;
            double maximum_relative_density_step = 0.0;
            std::string reference_branch;
            auto observe = [&](
                const std::vector<double>& composition,
                std::optional<detail::DensityRootIdentity>& previous_root
            ) {
                auto observation = endpoint.model.evaluate_fugacity_derivatives(
                    temperature_k,
                    pressure_pa,
                    composition,
                    "liquid",
                    parameters
                );
                maximum_root_residual = std::max(
                    maximum_root_residual,
                    std::abs(observation.root_residual_pa)
                );
                minimum_slope = std::min(
                    minimum_slope,
                    observation.pressure_density_derivative_pa_m3_per_mol
                );
                maximum_condition = std::max(
                    maximum_condition,
                    observation.density_condition_number
                );
                const double relative_bracket_width =
                    (observation.bracket_upper_mol_per_m3
                     - observation.bracket_lower_mol_per_m3)
                    / std::max(
                        1.0,
                        0.5
                            * (observation.bracket_upper_mol_per_m3
                               + observation.bracket_lower_mol_per_m3)
                    );
                maximum_relative_bracket_width = std::max(
                    maximum_relative_bracket_width,
                    relative_bracket_width
                );
                if (previous_root.has_value()) {
                    const double relative_density_step = std::abs(
                        observation.molar_density_mol_per_m3
                        - previous_root->molar_density
                    ) / std::max(1.0, previous_root->molar_density);
                    maximum_relative_density_step = std::max(
                        maximum_relative_density_step,
                        relative_density_step
                    );
                }
                const detail::DensityRootIdentity root{
                    observation.molar_density_mol_per_m3,
                    observation.stable_root_count,
                    observation.selected_stable_root_index,
                    observation.branch,
                };
                if (previous_root.has_value()
                    && !detail::same_continuous_density_root(
                        *previous_root, root, 1.0e-3
                    )) {
                    throw detail::DensityDomainError(
                        "neutral-reference density root is discontinuous"
                    );
                }
                previous_root = root;
                reference_branch = observation.branch;
                return observation;
            };
            std::optional<detail::DensityRootIdentity> reference_root;
            const auto reference_state = observe(
                reference_composition, reference_root
            );
            std::size_t row = 0;
            for (std::size_t component : reference_basis.neutral_components) {
                contractions[row] =
                    reference_state.log_fugacity_coefficient[component];
                pressure_derivatives[row] =
                    reference_state.pressure_derivative_per_pa[component];
                for (std::size_t parameter = 0;
                     parameter < request_count;
                     ++parameter) {
                    parameter_derivatives[
                        row * request_count + parameter
                    ] = reference_state.parameter_derivatives[
                        component * request_count + parameter
                    ];
                }
                ++row;
            }
            double reference_molality = nan;
            double reference_convergence_error = 0.0;
            double derivative_convergence_error = 0.0;
            auto append_pair = [&](std::size_t cation, std::size_t anion) {
                std::optional<detail::DensityRootIdentity> sequence_root =
                    reference_root;
                const int divisor = std::gcd(
                    std::abs(endpoint.component_charges[cation]),
                    std::abs(endpoint.component_charges[anion])
                );
                const std::array<int, 2> stoichiometry{
                    std::abs(endpoint.component_charges[anion]) / divisor,
                    std::abs(endpoint.component_charges[cation]) / divisor,
                };
                const auto sequence = detail::electrolyte_reference_sequence(
                    reference_composition,
                    endpoint.model.molar_masses_kg_per_mol(),
                    std::array<std::size_t, 2>{cation, anion},
                    stoichiometry,
                    [&](const std::vector<double>& composition) {
                        const auto observation = observe(
                            composition, sequence_root
                        );
                        std::vector<double> derivatives(
                            2 * (1 + request_count)
                        );
                        for (std::size_t ion = 0; ion < 2; ++ion) {
                            const std::size_t component =
                                ion == 0 ? cation : anion;
                            derivatives[ion * (1 + request_count)] =
                                observation.pressure_derivative_per_pa[
                                    component
                                ];
                            for (std::size_t parameter = 0;
                                 parameter < request_count;
                                 ++parameter) {
                                derivatives[
                                    ion * (1 + request_count) + 1 + parameter
                                ] = observation.parameter_derivatives[
                                    component * request_count + parameter
                                ];
                            }
                        }
                        return detail::ReferenceObservation{
                            {
                                observation.log_fugacity_coefficient[cation],
                                observation.log_fugacity_coefficient[anion],
                            },
                            std::move(derivatives),
                        };
                    }
                );
                if (!sequence.derivative.has_value()
                    || !sequence.derivative_convergence_error.has_value()) {
                    throw std::runtime_error(
                        "neutral-reference derivative limit is unavailable"
                    );
                }
                contractions[row] =
                    stoichiometry[0] * sequence.log_fugacity[0]
                    + stoichiometry[1] * sequence.log_fugacity[1];
                pressure_derivatives[row] =
                    stoichiometry[0] * (*sequence.derivative)[0]
                    + stoichiometry[1]
                        * (*sequence.derivative)[1 + request_count];
                for (std::size_t parameter = 0;
                     parameter < request_count;
                     ++parameter) {
                    parameter_derivatives[
                        row * request_count + parameter
                    ] = stoichiometry[0]
                            * (*sequence.derivative)[1 + parameter]
                        + stoichiometry[1]
                            * (*sequence.derivative)[
                                2 + request_count + parameter
                            ];
                }
                reference_molality = sequence.molality;
                reference_convergence_error = std::max(
                    reference_convergence_error, sequence.convergence_error
                );
                derivative_convergence_error = std::max(
                    derivative_convergence_error,
                    *sequence.derivative_convergence_error
                );
                ++row;
            };
            for (std::size_t anion : reference_basis.anions) {
                append_pair(reference_basis.cations.front(), anion);
            }
            for (std::size_t index = 1;
                 index < reference_basis.cations.size();
                 ++index) {
                append_pair(
                    reference_basis.cations[index],
                    reference_basis.anions.front()
                );
            }
            if (row != row_count || !std::isfinite(minimum_slope)
                || minimum_slope <= 0.0 || !std::isfinite(maximum_condition)
                || maximum_condition > 1.0e10
                || !std::isfinite(maximum_relative_bracket_width)
                || maximum_relative_bracket_width > 1.0e-8
                || !std::isfinite(maximum_relative_density_step)
                || !reference_root.has_value()
                || reference_branch.empty()) {
                throw detail::DensityDomainError(
                    "neutral-reference derivative path is ill-conditioned"
                );
            }
            const auto finite_values = [](const std::vector<double>& values) {
                return std::all_of(
                    values.begin(), values.end(), [](double item) {
                        return std::isfinite(item);
                    }
                );
            };
            if (!finite_values(basis) || !finite_values(contractions)
                || !finite_values(pressure_derivatives)
                || !finite_values(parameter_derivatives)
                || !std::isfinite(reference_molality)
                || !std::isfinite(reference_convergence_error)
                || !std::isfinite(derivative_convergence_error)) {
                throw std::runtime_error(
                    "neutral-reference derivative result is nonfinite"
                );
            }
            std::copy(basis.begin(), basis.end(), value.neutral_basis);
            std::copy(
                contractions.begin(),
                contractions.end(),
                value.log_fugacity_contractions
            );
            std::copy(
                reference_composition.begin(),
                reference_composition.end(),
                value.reference_composition
            );
            std::copy(
                pressure_derivatives.begin(),
                pressure_derivatives.end(),
                result->pressure_derivatives_per_pa
            );
            std::copy(
                parameter_derivatives.begin(),
                parameter_derivatives.end(),
                result->parameter_derivatives
            );
            value.derivative_availability =
                EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_NONE_V1;
            value.temperature_k = temperature_k;
            value.pressure_pa = pressure_pa;
            value.solvent_molar_mass_kg_per_mol = std::inner_product(
                reference_composition.begin(),
                reference_composition.end(),
                endpoint.model.molar_masses_kg_per_mol().begin(),
                0.0
            );
            value.reference_amount_mol = 1.0;
            value.reference_number_density_mol_per_m3 = 1.0;
            value.reference_molality_mol_per_kg = reference_molality;
            value.reference_convergence_error = reference_convergence_error;
            value.error[0] = '\0';
            std::snprintf(
                value.parameter_fingerprint,
                sizeof(value.parameter_fingerprint),
                "%s",
                endpoint.parameter_fingerprint.c_str()
            );
            std::snprintf(
                value.helmholtz_basis_id,
                sizeof(value.helmholtz_basis_id),
                "%s",
                EPCSAFT_NATIVE_HELMHOLTZ_BASIS_ID_V1
            );
            value.status = EPCSAFT_NATIVE_STATUS_OK_V1;
            result->derivative_availability =
                EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_PRESSURE_V1
                | (request_count == 0
                       ? EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_NONE_V1
                       : EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_PARAMETERS_V1);
            result->maximum_root_residual_pa = maximum_root_residual;
            result->source_pressure_min_pa = endpoint.pressure_min_pa;
            result->source_pressure_max_pa = endpoint.pressure_max_pa;
            result->minimum_pressure_density_derivative_pa_m3_per_mol =
                minimum_slope;
            result->maximum_density_condition_number = maximum_condition;
            result->reference_derivative_convergence_error =
                derivative_convergence_error;
            result->maximum_relative_root_bracket_width =
                maximum_relative_bracket_width;
            result->maximum_relative_root_density_step =
                maximum_relative_density_step;
            result->stable_root_count = reference_root->stable_root_count;
            result->selected_stable_root_index =
                reference_root->selected_stable_root_index;
            std::snprintf(
                result->reference_branch,
                sizeof(result->reference_branch),
                "%s",
                reference_branch.c_str()
            );
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
                "unknown neutral-reference derivative failure"
            );
        }
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_reacting_phase_parameters(
        const void* context,
        const char* expected_fingerprint,
        const char* expected_topology_fingerprint,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double volume_m3,
        double packing_fraction_min,
        double packing_fraction_max,
        const epcsaft_active_parameter_request_v1* requests,
        size_t request_count,
        epcsaft_reacting_phase_parameter_result_v1* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if (result->struct_size
            != sizeof(epcsaft_reacting_phase_parameter_result_v1)) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        const double nan = std::numeric_limits<double>::quiet_NaN();
        result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
        result->helmholtz_over_rt_reference_amount = nan;
        result->pressure_pa = nan;
        result->packing_fraction = nan;
        std::fill(
            std::begin(result->molar_volume_bounds_m3_per_mol),
            std::end(result->molar_volume_bounds_m3_per_mol),
            nan
        );
        result->parameter_fingerprint[0] = '\0';
        result->topology_fingerprint[0] = '\0';
        result->helmholtz_basis_id[0] = '\0';
        result->error[0] = '\0';
        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr || expected_fingerprint == nullptr
            || expected_topology_fingerprint == nullptr
            || amounts_mol == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "reacting phase required pointer is null"
            );
        }
        const auto& endpoint =
            *static_cast<const NativeSdkEndpoint*>(context);
        if (!endpoint.phase_callbacks_admitted
            || !endpoint.model.supports_electrolyte_native_sdk_v1()
            || endpoint.topology_fingerprint.empty()
            || std::strcmp(
                expected_fingerprint,
                endpoint.parameter_fingerprint.c_str()
            ) != 0
            || std::strcmp(
                expected_topology_fingerprint,
                endpoint.topology_fingerprint.c_str()
            ) != 0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1,
                "reacting phase model identity is unsupported"
            );
        }
        if (!std::isfinite(temperature_k)
            || temperature_k < endpoint.temperature_min_k
            || temperature_k > endpoint.temperature_max_k) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1,
                "reacting phase temperature is outside the source domain"
            );
        }
        const std::size_t state_count = endpoint.component_ids.size() + 1;
        if (amount_count != endpoint.component_ids.size()
            || result->state_coordinate_count != state_count
            || result->active_parameter_count != request_count
            || result->state_gradient_capacity != state_count
            || result->state_hessian_capacity != state_count * state_count
            || result->parameter_gradient_capacity != request_count
            || result->state_parameter_capacity != state_count * request_count
            || result->packing_parameter_capacity != request_count
            || result->packing_state_gradient_capacity != state_count
            || result->packing_state_hessian_capacity
                != state_count * state_count
            || result->molar_volume_bound_parameter_capacity
                != EPCSAFT_NATIVE_MOLAR_VOLUME_BOUND_COUNT_V1 * request_count
            || result->pressure_parameter_capacity != request_count
            || result->chemical_potential_parameter_capacity
                != endpoint.component_ids.size() * request_count
            || result->state_gradient == nullptr
            || result->state_hessian == nullptr
            || result->parameter_gradient == nullptr
            || result->state_parameter_derivatives == nullptr
            || result->packing_parameter_derivatives == nullptr
            || result->packing_state_gradient == nullptr
            || result->packing_state_hessian == nullptr
            || result->molar_volume_bound_parameter_derivatives_m3_per_mol
                == nullptr
            || result->pressure_parameter_derivatives_pa == nullptr
            || result->chemical_potential_parameter_derivatives_over_rt
                == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1,
                "reacting phase buffers do not match immutable dimensions"
            );
        }
        try {
            const auto parameters = active_parameters(
                endpoint, requests, request_count
            );
            const std::vector<double> amounts(
                amounts_mol, amounts_mol + amount_count
            );
            double total_amount = 0.0;
            for (const double amount : amounts) {
                total_amount += amount;
            }
            std::vector<double> composition(amount_count);
            for (std::size_t component = 0;
                 component < amount_count;
                 ++component) {
                composition[component] = amounts[component] / total_amount;
            }
            const auto packing =
                endpoint.model.evaluate_packing_fraction_parameters(
                    temperature_k,
                    amounts,
                    volume_m3,
                    endpoint.ion_mole_fraction_max,
                    parameters
                );
            const auto phase = endpoint.model.evaluate_reacting_phase_parameters(
                temperature_k,
                amounts,
                volume_m3,
                endpoint.ion_mole_fraction_max,
                parameters
            );
            const auto volume_bounds =
                endpoint.model.molar_volume_bounds_parameters(
                    temperature_k,
                    composition,
                    packing_fraction_min,
                    packing_fraction_max,
                    endpoint.ion_mole_fraction_max,
                    parameters
                );
            const std::size_t joint_count = state_count + request_count;
            if (phase.gradient.size() != joint_count
                || phase.hessian.size() != joint_count * joint_count
                || packing.gradient.size() != joint_count
                || packing.hessian.size() != joint_count * joint_count) {
                throw std::runtime_error(
                    "reacting phase tensor dimensions changed"
                );
            }
            std::vector<double> state_gradient(state_count);
            std::vector<double> state_hessian(state_count * state_count);
            std::vector<double> parameter_gradient(request_count);
            std::vector<double> mixed(state_count * request_count);
            std::vector<double> packing_parameter(request_count);
            std::vector<double> packing_state_gradient(state_count);
            std::vector<double> packing_state_hessian(
                state_count * state_count
            );
            std::vector<double> pressure_parameter(request_count);
            std::vector<double> chemical_parameter(
                amount_count * request_count
            );
            for (std::size_t state = 0; state < state_count; ++state) {
                state_gradient[state] = phase.gradient[state];
                packing_state_gradient[state] = packing.gradient[state];
                for (std::size_t other = 0; other < state_count; ++other) {
                    state_hessian[state * state_count + other] =
                        phase.hessian[state * joint_count + other];
                    packing_state_hessian[state * state_count + other] =
                        packing.hessian[state * joint_count + other];
                }
                for (std::size_t parameter = 0;
                     parameter < request_count;
                     ++parameter) {
                    mixed[state * request_count + parameter] =
                        phase.hessian[
                            state * joint_count + state_count + parameter
                        ];
                }
            }
            for (std::size_t parameter = 0;
                 parameter < request_count;
                 ++parameter) {
                parameter_gradient[parameter] =
                    phase.gradient[state_count + parameter];
                packing_parameter[parameter] =
                    packing.gradient[state_count + parameter];
                pressure_parameter[parameter] =
                    -detail::GAS_CONSTANT * temperature_k
                    * mixed[
                        (state_count - 1) * request_count + parameter
                    ];
                for (std::size_t component = 0;
                     component < amount_count;
                     ++component) {
                    chemical_parameter[
                        component * request_count + parameter
                    ] = mixed[component * request_count + parameter];
                }
            }
            const auto finite_values = [](const std::vector<double>& values) {
                return std::all_of(
                    values.begin(), values.end(), [](double value) {
                        return std::isfinite(value);
                    }
                );
            };
            if (!finite_values(state_gradient)
                || !finite_values(state_hessian)
                || !finite_values(parameter_gradient)
                || !finite_values(mixed)
                || !finite_values(packing_parameter)
                || !finite_values(packing_state_gradient)
                || !finite_values(packing_state_hessian)
                || !finite_values(
                    volume_bounds.parameter_derivatives_m3_per_mol
                )
                || !finite_values(pressure_parameter)
                || !finite_values(chemical_parameter)
                || !std::isfinite(phase.helmholtz_over_rt_reference_amount)
                || !std::isfinite(phase.pressure_pa)
                || !std::isfinite(packing.packing_fraction)) {
                throw std::runtime_error(
                    "reacting phase result contains nonfinite values"
                );
            }
            std::copy(
                state_gradient.begin(),
                state_gradient.end(),
                result->state_gradient
            );
            std::copy(
                state_hessian.begin(),
                state_hessian.end(),
                result->state_hessian
            );
            std::copy(
                parameter_gradient.begin(),
                parameter_gradient.end(),
                result->parameter_gradient
            );
            std::copy(
                mixed.begin(),
                mixed.end(),
                result->state_parameter_derivatives
            );
            std::copy(
                packing_parameter.begin(),
                packing_parameter.end(),
                result->packing_parameter_derivatives
            );
            std::copy(
                packing_state_gradient.begin(),
                packing_state_gradient.end(),
                result->packing_state_gradient
            );
            std::copy(
                packing_state_hessian.begin(),
                packing_state_hessian.end(),
                result->packing_state_hessian
            );
            std::copy(
                volume_bounds.parameter_derivatives_m3_per_mol.begin(),
                volume_bounds.parameter_derivatives_m3_per_mol.end(),
                result
                    ->molar_volume_bound_parameter_derivatives_m3_per_mol
            );
            std::copy(
                pressure_parameter.begin(),
                pressure_parameter.end(),
                result->pressure_parameter_derivatives_pa
            );
            std::copy(
                chemical_parameter.begin(),
                chemical_parameter.end(),
                result->chemical_potential_parameter_derivatives_over_rt
            );
            result->helmholtz_over_rt_reference_amount =
                phase.helmholtz_over_rt_reference_amount;
            result->pressure_pa = phase.pressure_pa;
            result->packing_fraction = packing.packing_fraction;
            std::copy(
                volume_bounds.bounds_m3_per_mol.begin(),
                volume_bounds.bounds_m3_per_mol.end(),
                result->molar_volume_bounds_m3_per_mol
            );
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
                "unknown reacting phase parameter failure"
            );
        }
    }

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
        const auto& reference_basis = endpoint.neutral_reference_basis;
        const std::size_t row_count = reference_basis.row_count();
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
            if (row_count != component_count - 1u) {
                return fail(
                    EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1,
                    "neutral reference component topology is incomplete"
                );
            }
            const std::vector<double> reference_composition =
                reference_basis.reference_composition;
            const auto reference_state = endpoint.model.evaluate_pressure(
                temperature_k, pressure_pa, reference_composition, "liquid"
            );
            if (!reference_state.log_fugacity_coefficient.has_value()) {
                throw std::runtime_error("neutral reference fugacity is unavailable");
            }
            std::vector<double> basis = reference_basis.rows;
            std::vector<double> contractions(row_count, 0.0);
            std::size_t row = 0;
            for (std::size_t component : reference_basis.neutral_components) {
                contractions[row] =
                    (*reference_state.log_fugacity_coefficient)[component];
                ++row;
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
            for (std::size_t anion : reference_basis.anions) {
                append_pair(reference_basis.cations.front(), anion);
            }
            for (std::size_t index = 1;
                 index < reference_basis.cations.size();
                 ++index) {
                append_pair(
                    reference_basis.cations[index],
                    reference_basis.anions.front()
                );
            }
            if (row != row_count) {
                throw std::runtime_error("neutral reference basis dimension changed");
            }
            const auto& molar_masses = endpoint.model.molar_masses_kg_per_mol();
            const double solvent_molar_mass = std::inner_product(
                reference_composition.begin(),
                reference_composition.end(),
                molar_masses.begin(),
                0.0
            );
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
    static int evaluate_fixed_topology_association(
        const void* context,
        const char* expected_parameter_fingerprint,
        const char* expected_topology_fingerprint,
        const char* expected_artifact_fingerprint,
        double temperature_k,
        double amount_mol,
        double volume_m3,
        const epcsaft_fixed_topology_parameter_request_v1* requests,
        size_t request_count,
        epcsaft_fixed_topology_phase_result_v1* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if (result->struct_size
            != sizeof(epcsaft_fixed_topology_phase_result_v1)) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        result->status = EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        result->error[0] = '\0';
        result->parameter_fingerprint[0] = '\0';
        result->topology_fingerprint[0] = '\0';
        result->artifact_fingerprint[0] = '\0';
        result->helmholtz_basis_id[0] = '\0';
        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr || expected_parameter_fingerprint == nullptr
            || expected_topology_fingerprint == nullptr
            || expected_artifact_fingerprint == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "model context and expected fingerprints are required"
            );
        }
        const auto& endpoint = *static_cast<const NativeSdkEndpoint*>(context);
        if (endpoint.fixed_topology_association.struct_size == 0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1,
                "model does not advertise fixed-topology association parameters"
            );
        }
        if (std::strcmp(
                expected_parameter_fingerprint,
                endpoint.parameter_fingerprint.c_str()
            ) != 0
            || std::strcmp(
                expected_topology_fingerprint,
                endpoint.topology_fingerprint.c_str()
            ) != 0
            || std::strcmp(
                expected_artifact_fingerprint,
                endpoint.artifact_fingerprint.c_str()
            ) != 0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1,
                "fixed-topology association fingerprint mismatch"
            );
        }
        if (request_count == 0 || requests == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "a nonempty active parameter request is required"
            );
        }
        if (request_count
            > endpoint.fixed_topology_parameter_slots.size()) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "active parameter request exceeds the advertised slot count"
            );
        }
        const std::size_t coordinate_count = 2 + request_count;
        if (result->coordinate_count != coordinate_count
            || result->active_parameter_count != request_count
            || result->gradient_capacity != coordinate_count
            || result->hessian_capacity != coordinate_count * coordinate_count
            || result->gradient == nullptr || result->hessian == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1,
                "fixed-topology result dimensions or buffers do not match"
            );
        }
        if (!std::isfinite(temperature_k) || temperature_k <= 0.0
            || !std::isfinite(amount_mol) || amount_mol <= 0.0
            || !std::isfinite(volume_m3) || volume_m3 <= 0.0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "fixed-topology phase state must be positive and finite"
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
            std::vector<NativeFixedTopologyParameter> active_parameters;
            active_parameters.reserve(request_count);
            std::size_t previous_slot_index = 0;
            bool has_previous_slot = false;
            for (std::size_t index = 0; index < request_count; ++index) {
                const auto& request = requests[index];
                if (request.struct_size
                        != sizeof(epcsaft_fixed_topology_parameter_request_v1)
                    || request.unit == nullptr
                    || !std::isfinite(request.value)) {
                    throw std::invalid_argument(
                        "fixed-topology active parameter request is malformed"
                    );
                }
                if (request.value <= 0.0) {
                    return fail(
                        EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1,
                        "fixed-topology active parameter value is outside its "
                        "advertised physical domain"
                    );
                }
                const auto slot = std::find_if(
                    endpoint.fixed_topology_parameter_slots.begin(),
                    endpoint.fixed_topology_parameter_slots.end(),
                    [&](const auto& candidate) {
                        return candidate.parameter_family
                                == request.parameter_family
                            && candidate.component_index
                                == request.component_index
                            && candidate.site_index_a == request.site_index_a
                            && candidate.site_index_b == request.site_index_b
                            && std::strcmp(candidate.unit, request.unit) == 0;
                    }
                );
                if (slot == endpoint.fixed_topology_parameter_slots.end()) {
                    throw std::invalid_argument(
                        "active parameter is not advertised by the fixed topology"
                    );
                }
                const auto slot_index = static_cast<std::size_t>(std::distance(
                    endpoint.fixed_topology_parameter_slots.begin(), slot
                ));
                if (has_previous_slot
                    && slot_index <= previous_slot_index) {
                    throw std::invalid_argument(
                        "active parameters must follow advertised slot order"
                    );
                }
                previous_slot_index = slot_index;
                has_previous_slot = true;
                NativeFixedTopologyParameterKind kind;
                switch (request.parameter_family) {
                    case EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_COUNT_V1:
                        kind = NativeFixedTopologyParameterKind::SegmentCount;
                        break;
                    case EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_DIAMETER_V1:
                        kind = NativeFixedTopologyParameterKind::SegmentDiameter;
                        break;
                    case EPCSAFT_NATIVE_PARAMETER_FAMILY_DISPERSION_ENERGY_OVER_K_V1:
                        kind =
                            NativeFixedTopologyParameterKind::DispersionEnergyOverK;
                        break;
                    case EPCSAFT_NATIVE_PARAMETER_FAMILY_ASSOCIATION_ENERGY_OVER_K_V1:
                        kind =
                            NativeFixedTopologyParameterKind::AssociationEnergyOverK;
                        break;
                    case EPCSAFT_NATIVE_PARAMETER_FAMILY_ASSOCIATION_VOLUME_V1:
                        kind =
                            NativeFixedTopologyParameterKind::AssociationVolume;
                        break;
                    default:
                        throw std::logic_error(
                            "advertised fixed-topology family is unsupported"
                        );
                }
                active_parameters.push_back({
                    kind,
                    request.site_index_a < 0
                        ? 0u
                        : static_cast<std::size_t>(request.site_index_a),
                    request.site_index_b < 0
                        ? 0u
                        : static_cast<std::size_t>(request.site_index_b),
                    request.value,
                });
            }
            const auto block = endpoint.model.evaluate_fixed_topology_association(
                temperature_k,
                amount_mol,
                volume_m3,
                active_parameters
            );
            if (block.gradient.size() != coordinate_count
                || block.hessian.size() != coordinate_count * coordinate_count
                || !std::isfinite(
                    block.helmholtz_over_rt_reference_amount
                )
                || !std::isfinite(block.pressure_pa)
                || !std::isfinite(block.chemical_potential_over_rt)
                || !std::all_of(
                    block.gradient.begin(), block.gradient.end(),
                    [](double value) { return std::isfinite(value); }
                )
                || !std::all_of(
                    block.hessian.begin(), block.hessian.end(),
                    [](double value) { return std::isfinite(value); }
                )) {
                throw std::runtime_error(
                    "fixed-topology phase result is incomplete or nonfinite"
                );
            }
            std::copy(
                block.gradient.begin(), block.gradient.end(), result->gradient
            );
            std::copy(
                block.hessian.begin(), block.hessian.end(), result->hessian
            );
            result->helmholtz_over_rt_reference_amount =
                block.helmholtz_over_rt_reference_amount;
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
            std::snprintf(
                result->artifact_fingerprint,
                sizeof(result->artifact_fingerprint),
                "%s",
                endpoint.artifact_fingerprint.c_str()
            );
            std::snprintf(
                result->helmholtz_basis_id,
                sizeof(result->helmholtz_basis_id),
                "%s",
                EPCSAFT_NATIVE_SDK_V1_HELMHOLTZ_BASIS_ID
            );
            result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
            return result->status;
        } catch (const detail::PolarDomainError& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1, error.what());
        } catch (const std::logic_error& error) {
            return fail(EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1, error.what());
        } catch (const std::exception& error) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1, error.what());
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1,
                "unknown fixed-topology association evaluation failure"
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
    static int evaluate_electrolyte_value(
        const void* context,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double volume_m3,
        epcsaft_electrolyte_phase_value_result_v1* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if (result->struct_size
            != sizeof(epcsaft_electrolyte_phase_value_result_v1)) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        const uint32_t result_size = result->struct_size;
        std::memset(result, 0, sizeof(*result));
        result->struct_size = result_size;
        result->helmholtz_over_rt_reference_amount =
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
        const auto& endpoint = *static_cast<const NativeSdkEndpoint*>(context);
        if (!endpoint.phase_callbacks_admitted) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1,
                "model is not admitted for this phase callback"
            );
        }
        if (amounts_mol == nullptr
            || amount_count != endpoint.table.component_count
            || !std::isfinite(temperature_k) || temperature_k <= 0.0
            || !std::isfinite(volume_m3) || volume_m3 <= 0.0
            || !std::all_of(
                amounts_mol,
                amounts_mol + amount_count,
                [](double amount) {
                    return std::isfinite(amount) && amount > 0.0;
                }
            )) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "temperature, amounts, and volume must be positive and finite"
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
            const std::vector<double> amounts(
                amounts_mol, amounts_mol + amount_count
            );
            result->helmholtz_over_rt_reference_amount =
                endpoint.model.evaluate_electrolyte_phase_value(
                    temperature_k,
                    amounts,
                    volume_m3,
                    endpoint.ion_mole_fraction_max
                );
            std::snprintf(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                "%s",
                endpoint.parameter_fingerprint.c_str()
            );
            return result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
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
    const std::string& topology_fingerprint,
    const std::vector<std::string>& site_ids,
    const std::string& artifact_fingerprint
) {
    const bool fixed_topology_admitted =
        fixed_topology_site_ids_valid(model, site_ids)
        && sha256_fingerprint_valid(artifact_fingerprint);
    if (!active_born_admitted && !aqueous_kij_admitted
        && !ion_fraction_suppression_admitted
        && !ion_solvation_kij_admitted
        && !neutral_reference_admitted && !aqueous_solvation_factor_admitted
        && (!model.supports_native_sdk_v1()
        && !fixed_topology_admitted
        && !model.supports_mixture_native_sdk_v1()
        && !model.supports_electrolyte_native_sdk_v1())) {
        throw std::invalid_argument(
            "native SDK v1 requires a supported phase-block model"
        );
    }
    if (active_born_admitted
        && !active_born_model_admitted(
            component_ids,
            model.charge_numbers()
        )) {
        throw std::invalid_argument("unsupported active Born model identity");
    }
    if (aqueous_kij_admitted
        && !aqueous_kij_model_admitted(
            component_ids,
            model.charge_numbers()
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
            component_ids, model.charge_numbers()
        )) {
        throw std::invalid_argument(
            "unsupported active ion-solvation k_ij model identity"
        );
    }
    if (aqueous_solvation_factor_admitted
        && !aqueous_solvation_factor_model_admitted(
            component_ids,
            model.charge_numbers()
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
    if (model.supports_fixed_topology_association_native_sdk_v1()
        && !site_ids.empty() && !fixed_topology_admitted) {
        throw std::invalid_argument(
            "invalid association site or installed-artifact identity"
        );
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
        site_ids,
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
        topology_fingerprint,
        artifact_fingerprint
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
