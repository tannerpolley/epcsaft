#include "owner_evaluator.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace epcsaft_engine::equilibrium {
namespace {

constexpr const char* kHomogeneousKind =
    "D-034.fixed_TP_homogeneous_reactive";
constexpr const char* kBubbleVleKind =
    "D-035.fixed_topology_reactive_bubble_vle";

interfaces::EvaluationFailure failure(
    std::string code,
    std::string reason
) {
    return {std::move(code), std::move(reason)};
}

bool finite(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
    });
}

bool sha256_fingerprint(const std::string& value) {
    return value.size() == 71 && value.rfind("sha256:", 0) == 0
        && std::all_of(value.begin() + 7, value.end(),
            [](unsigned char character) { return std::isxdigit(character) != 0; });
}

std::string sha256_digest(const std::string& input) {
    constexpr std::array<std::uint32_t, 64> constants{
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
        0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
        0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
        0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
        0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
        0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
        0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
        0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
    };
    std::array<std::uint32_t, 8> state{
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
        0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
    };
    std::vector<unsigned char> bytes(input.begin(), input.end());
    const std::uint64_t bit_count = static_cast<std::uint64_t>(bytes.size()) * 8;
    bytes.push_back(0x80);
    while ((bytes.size() + 8) % 64 != 0) {
        bytes.push_back(0);
    }
    for (int shift = 56; shift >= 0; shift -= 8) {
        bytes.push_back(static_cast<unsigned char>(bit_count >> shift));
    }
    const auto rotate_right = [](std::uint32_t value, unsigned amount) {
        return (value >> amount) | (value << (32 - amount));
    };
    for (std::size_t offset = 0; offset < bytes.size(); offset += 64) {
        std::array<std::uint32_t, 64> words{};
        for (std::size_t index = 0; index < 16; ++index) {
            const std::size_t position = offset + index * 4;
            words[index] = (static_cast<std::uint32_t>(bytes[position]) << 24)
                | (static_cast<std::uint32_t>(bytes[position + 1]) << 16)
                | (static_cast<std::uint32_t>(bytes[position + 2]) << 8)
                | static_cast<std::uint32_t>(bytes[position + 3]);
        }
        for (std::size_t index = 16; index < words.size(); ++index) {
            const auto s0 = rotate_right(words[index - 15], 7)
                ^ rotate_right(words[index - 15], 18)
                ^ (words[index - 15] >> 3);
            const auto s1 = rotate_right(words[index - 2], 17)
                ^ rotate_right(words[index - 2], 19)
                ^ (words[index - 2] >> 10);
            words[index] = words[index - 16] + s0 + words[index - 7] + s1;
        }
        std::uint32_t a = state[0];
        std::uint32_t b = state[1];
        std::uint32_t c = state[2];
        std::uint32_t d = state[3];
        std::uint32_t e = state[4];
        std::uint32_t f = state[5];
        std::uint32_t g = state[6];
        std::uint32_t h = state[7];
        for (std::size_t index = 0; index < words.size(); ++index) {
            const auto sum1 = rotate_right(e, 6) ^ rotate_right(e, 11)
                ^ rotate_right(e, 25);
            const auto choice = (e & f) ^ ((~e) & g);
            const auto temporary1 = h + sum1 + choice + constants[index]
                + words[index];
            const auto sum0 = rotate_right(a, 2) ^ rotate_right(a, 13)
                ^ rotate_right(a, 22);
            const auto majority = (a & b) ^ (a & c) ^ (b & c);
            const auto temporary2 = sum0 + majority;
            h = g;
            g = f;
            f = e;
            e = d + temporary1;
            d = c;
            c = b;
            b = a;
            a = temporary1 + temporary2;
        }
        state[0] += a;
        state[1] += b;
        state[2] += c;
        state[3] += d;
        state[4] += e;
        state[5] += f;
        state[6] += g;
        state[7] += h;
    }
    std::ostringstream result;
    result << std::hex << std::setfill('0');
    for (const auto word : state) {
        result << std::setw(8) << word;
    }
    return result.str();
}

void append_field(std::string& payload, const std::string& value) {
    payload += std::to_string(value.size());
    payload.push_back(':');
    payload += value;
    payload.push_back('|');
}

void append_json_string(std::string& payload, const std::string& value) {
    constexpr char digits[] = "0123456789abcdef";
    const auto append_hex4 = [&](std::uint16_t code_unit) {
        payload += "\\u";
        payload.push_back(digits[(code_unit >> 12) & 0x0f]);
        payload.push_back(digits[(code_unit >> 8) & 0x0f]);
        payload.push_back(digits[(code_unit >> 4) & 0x0f]);
        payload.push_back(digits[code_unit & 0x0f]);
    };
    payload.push_back('"');
    for (std::size_t index = 0; index < value.size();) {
        const unsigned char character = static_cast<unsigned char>(value[index]);
        if (character >= 0x80) {
            std::uint32_t codepoint = 0;
            std::size_t width = 0;
            if (character >= 0xc2 && character <= 0xdf) {
                codepoint = character & 0x1f;
                width = 2;
            } else if (character >= 0xe0 && character <= 0xef) {
                codepoint = character & 0x0f;
                width = 3;
            } else if (character >= 0xf0 && character <= 0xf4) {
                codepoint = character & 0x07;
                width = 4;
            } else {
                throw std::invalid_argument("source specification string is not valid UTF-8");
            }
            if (index + width > value.size()) {
                throw std::invalid_argument("source specification string is not valid UTF-8");
            }
            for (std::size_t continuation = 1; continuation < width; ++continuation) {
                const unsigned char next = static_cast<unsigned char>(
                    value[index + continuation]
                );
                if ((next & 0xc0) != 0x80) {
                    throw std::invalid_argument(
                        "source specification string is not valid UTF-8"
                    );
                }
                codepoint = (codepoint << 6) | (next & 0x3f);
            }
            if ((width == 3 && codepoint < 0x800)
                || (width == 4 && codepoint < 0x10000)
                || codepoint > 0x10ffff
                || (codepoint >= 0xd800 && codepoint <= 0xdfff)) {
                throw std::invalid_argument("source specification string is not valid UTF-8");
            }
            if (codepoint <= 0xffff) {
                append_hex4(static_cast<std::uint16_t>(codepoint));
            } else {
                codepoint -= 0x10000;
                append_hex4(static_cast<std::uint16_t>(0xd800 + (codepoint >> 10)));
                append_hex4(static_cast<std::uint16_t>(0xdc00 + (codepoint & 0x3ff)));
            }
            index += width;
            continue;
        }
        switch (character) {
        case '"': payload += "\\\""; break;
        case '\\': payload += "\\\\"; break;
        case '\b': payload += "\\b"; break;
        case '\f': payload += "\\f"; break;
        case '\n': payload += "\\n"; break;
        case '\r': payload += "\\r"; break;
        case '\t': payload += "\\t"; break;
        default:
            if (character < 0x20) {
                payload += "\\u00";
                payload.push_back(digits[character >> 4]);
                payload.push_back(digits[character & 0x0f]);
            } else {
                payload.push_back(static_cast<char>(character));
            }
        }
        ++index;
    }
    payload.push_back('"');
}

void append_json_number(std::string& payload, double value) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument("source specification number must be finite");
    }
    std::array<char, 64> buffer{};
    const auto [end, error] = std::to_chars(
        buffer.data(), buffer.data() + buffer.size(), value,
        std::chars_format::general
    );
    if (error != std::errc()) {
        throw std::invalid_argument("source specification number is not representable");
    }
    std::string rendered(buffer.data(), end);
    const std::size_t exponent_marker = rendered.find_first_of("eE");
    if (exponent_marker != std::string::npos) {
        std::string mantissa = rendered.substr(0, exponent_marker);
        bool exponent_negative = false;
        std::size_t exponent_start = exponent_marker + 1;
        if (rendered[exponent_start] == '+' || rendered[exponent_start] == '-') {
            exponent_negative = rendered[exponent_start] == '-';
            ++exponent_start;
        }
        int exponent = 0;
        for (; exponent_start < rendered.size(); ++exponent_start) {
            exponent = exponent * 10 + (rendered[exponent_start] - '0');
        }
        if (exponent_negative) {
            exponent = -exponent;
        }
        // Python's JSON encoder uses fixed notation for -4 <= exponent < 16.
        if (exponent >= -4 && exponent < 16) {
            bool negative = !mantissa.empty() && mantissa.front() == '-';
            if (negative) {
                mantissa.erase(0, 1);
            }
            const std::size_t dot = mantissa.find('.');
            const int decimal_position = static_cast<int>(
                dot == std::string::npos ? mantissa.size() : dot
            );
            if (dot != std::string::npos) {
                mantissa.erase(dot, 1);
            }
            const int target = decimal_position + exponent;
            if (target <= 0) {
                rendered = "0." + std::string(static_cast<std::size_t>(-target), '0')
                    + mantissa;
            } else if (target >= static_cast<int>(mantissa.size())) {
                rendered = mantissa + std::string(
                    static_cast<std::size_t>(target - mantissa.size()), '0'
                );
            } else {
                rendered = mantissa.substr(0, static_cast<std::size_t>(target))
                    + "." + mantissa.substr(static_cast<std::size_t>(target));
            }
            if (rendered.find('.') == std::string::npos) {
                rendered += ".0";
            }
            if (negative) {
                rendered.insert(0, 1, '-');
            }
        } else {
            std::string exponent_digits = std::to_string(
                exponent < 0 ? -exponent : exponent
            );
            if (exponent_digits.size() < 2) {
                exponent_digits.insert(0, 1, '0');
            }
            rendered = mantissa + "e" + (exponent < 0 ? "-" : "+")
                + exponent_digits;
        }
    } else if (rendered.find('.') == std::string::npos) {
        rendered += ".0";
    }
    payload += rendered;
}

template <typename Append>
void append_json_array(std::string& payload, std::size_t count, Append append) {
    payload.push_back('[');
    for (std::size_t index = 0; index < count; ++index) {
        if (index != 0) {
            payload.push_back(',');
        }
        append(index);
    }
    payload.push_back(']');
}

std::string canonical_source_specification_fingerprint(
    const GrepeRequest& request,
    bool active_parameters_requested
) {
    if (!request.source.has_value()) {
        return {};
    }
    const GrepeSourceStandardState& source = *request.source;
    const epcsaft::NativeSourceReferenceRequest& reference = source.reference;
    const std::uint32_t expected_derivatives =
        2u | (active_parameters_requested ? 8u : 0u);
    if (reference.required_derivatives != expected_derivatives) {
        throw std::invalid_argument("source specification derivative contract changed");
    }
    const auto component_identity = [&](std::size_t index) -> const std::string& {
        if (index >= reference.component_ids.size()) {
            throw std::invalid_argument("source specification component index is invalid");
        }
        return reference.component_ids[index];
    };
    std::string payload = "{";
    payload += "\"activity_convention_id\":";
    append_json_string(payload, reference.activity_convention_id);
    payload += ",\"component_ids\":";
    append_json_array(payload, reference.component_ids.size(), [&](std::size_t index) {
        append_json_string(payload, reference.component_ids[index]);
    });
    payload += ",\"ion_pairs\":";
    if (reference.ion_pair_indices.size() != reference.ion_pair_stoichiometry.size()) {
        throw std::invalid_argument("source specification ion-pair metadata is inconsistent");
    }
    append_json_array(payload, reference.ion_pair_indices.size(), [&](std::size_t index) {
        payload.push_back('[');
        append_json_array(payload, 2, [&](std::size_t member) {
            append_json_string(payload, component_identity(
                reference.ion_pair_indices[index][member]
            ));
        });
        payload.push_back(',');
        append_json_array(payload, 2, [&](std::size_t member) {
            payload += std::to_string(
                reference.ion_pair_stoichiometry[index][member]
            );
        });
        payload.push_back(']');
    });
    payload += ",\"phase\":";
    append_json_string(payload, reference.phase);
    payload += ",\"pure_components\":";
    if (reference.pure_component_indices.size() != reference.pure_component_phases.size()) {
        throw std::invalid_argument("source specification pure-component metadata is inconsistent");
    }
    append_json_array(payload, reference.pure_component_indices.size(), [&](std::size_t index) {
        payload.push_back('[');
        append_json_string(payload, component_identity(reference.pure_component_indices[index]));
        payload.push_back(',');
        append_json_string(payload, reference.pure_component_phases[index]);
        payload.push_back(']');
    });
    payload += ",\"reference_convention\":";
    append_json_string(payload, reference.reference_convention);
    payload += ",\"reference_pressure_pa\":";
    append_json_number(payload, reference.source_reference_pressure_pa);
    payload += ",\"reference_state_id\":";
    append_json_string(payload, reference.reference_state_id);
    payload += ",\"required_derivatives\":[";
    append_json_string(payload, "pressure");
    if (active_parameters_requested) {
        payload.push_back(',');
        append_json_string(payload, "parameters");
    }
    payload += "],\"schema\":\"epcsaft-source-reference-declaration-v1\",\"solvent_mole_fractions\":";
    append_json_array(payload, reference.source_solvent_mole_fractions.size(), [&](std::size_t index) {
        append_json_number(payload, reference.source_solvent_mole_fractions[index]);
    });
    payload += ",\"standard_molality_mol_per_kg\":";
    append_json_number(payload, reference.standard_molality_mol_per_kg);
    payload.push_back('}');
    return "sha256:" + sha256_digest(payload);
}

bool ends_with(const std::string& value, const char* suffix) {
    const std::size_t length = std::char_traits<char>::length(suffix);
    return value.size() >= length
        && value.compare(value.size() - length, length, suffix) == 0;
}

std::vector<std::string> canonical_parameter_units(
    const std::vector<epcsaft::NativeActiveCoordinate>& coordinates
) {
    std::vector<std::string> result;
    result.reserve(coordinates.size());
    for (const auto& coordinate : coordinates) {
        const auto& family = coordinate.family;
        if (family == "segment_diameter" || family == "born_diameter") {
            result.push_back("angstrom");
        } else if (family == "dispersion_energy_over_k"
                   || family == "association_energy_over_k"
                   || family == "k_ij_reciprocal_temperature_slope") {
            result.push_back("kelvin");
        } else if (family == "exponential_temperature_coefficient") {
            result.push_back("1 / kelvin");
        } else if (family.rfind("temperature_polynomial_coefficient_", 0) == 0) {
            result.push_back("1 / kelvin ** " + family.substr(35));
        } else if (family == "reaction_correlation") {
            result.push_back(ends_with(coordinate.identity, "_per_k") ? "1 / kelvin"
                : ends_with(coordinate.identity, "_k") ? "kelvin"
                : "dimensionless");
        } else if (family == "segment_count" || family == "relative_permittivity"
                   || family == "solvation_factor" || family == "association_volume"
                   || family == "k_ij" || family == "l_ij" || family == "k_hb_ij"
                   || family == "ion_fraction_suppression_coefficient"
                   || family == "ionic_region_relative_permittivity"
                   || family == "ion_specific_suppression_coefficient"
                   || family == "correlation_integral_parameter") {
            result.push_back("dimensionless");
        } else {
            throw std::invalid_argument("unsupported EOS parameter unit family");
        }
    }
    return result;
}

std::vector<std::string> canonical_parameter_bases(
    const std::vector<epcsaft::NativeActiveCoordinate>& coordinates
) {
    return std::vector<std::string>(coordinates.size(), "physical");
}

bool unique_nonempty(const std::vector<std::string>& values) {
    if (!std::all_of(values.begin(), values.end(), [](const std::string& value) {
        return !value.empty();
    })) {
        return false;
    }
    std::vector<std::string> sorted = values;
    std::sort(sorted.begin(), sorted.end());
    return std::adjacent_find(sorted.begin(), sorted.end()) == sorted.end();
}

std::string failure_code(const GrepeResult& result) {
    if (!result.failure_code.empty()) {
        return result.failure_code;
    }
    const auto& equilibrium = result.equilibrium;
    if (equilibrium.eos_domain_status == "failed") {
        return "eos_domain_rejection";
    }
    if (equilibrium.physical_status == "failed") {
        return "physical_rejection";
    }
    if (equilibrium.sensitivity_status == "unavailable") {
        return "singular_sensitivity";
    }
    return "numerical_convergence_failure";
}

interfaces::EvaluationFailure invalid_result(
    const std::exception& error
) {
    const std::string reason = error.what();
    const bool domain = reason.find("EOS") != std::string::npos
        || reason.find("packing") != std::string::npos
        || reason.find("domain") != std::string::npos;
    return failure(domain ? "eos_domain_rejection" : "owner_failure", reason);
}

epcsaft::NativeModel model_at(
    const epcsaft::NativeModel& model,
    const std::vector<double>& values
) {
    std::vector<epcsaft::NativeActiveCoordinate> coordinates = model.active_coordinates();
    if (coordinates.size() != values.size()) {
        throw std::invalid_argument("active parameter dimensions changed");
    }
    if (coordinates.empty()) {
        return model;
    }
    for (std::size_t index = 0; index < coordinates.size(); ++index) {
        coordinates[index].value = values[index];
    }
    return model.with_active_coordinates(coordinates);
}

std::vector<double> values_for_model(
    const epcsaft::NativeModel& model,
    const std::vector<EvaluatorParameter>& parameters,
    const std::vector<double>& values
) {
    std::vector<double> result;
    for (const auto& coordinate : model.active_coordinates()) {
        const auto found = std::find_if(
            parameters.begin(), parameters.end(), [&](const EvaluatorParameter& parameter) {
                return parameter.identity == coordinate.identity;
            }
        );
        if (found == parameters.end()) {
            throw std::invalid_argument("active parameter order changed");
        }
        result.push_back(values[static_cast<std::size_t>(found - parameters.begin())]);
    }
    return result;
}

std::vector<std::string> identities(
    const std::vector<epcsaft::NativeActiveCoordinate>& coordinates
) {
    return epcsaft::active_parameter_identities(coordinates);
}

std::string canonical_hash(const std::vector<std::string>& fields) {
    std::string payload = "epcsaft.equilibrium-owner.v1|";
    for (const auto& field : fields) {
        append_field(payload, field);
    }
    return "sha256:" + sha256_digest(payload);
}

std::string canonical_double(double value) {
    std::ostringstream rendered;
    rendered << std::setprecision(17) << value;
    return rendered.str();
}

void append_doubles(std::vector<std::string>& fields, const std::vector<double>& values) {
    for (const double value : values) {
        fields.push_back(canonical_double(value));
    }
}

void append_indices(std::vector<std::string>& fields, const std::vector<std::size_t>& values) {
    for (const std::size_t value : values) {
        fields.push_back(std::to_string(value));
    }
}

void append_matrix(std::vector<std::string>& fields, const DenseMatrix& matrix) {
    fields.push_back(std::to_string(matrix.rows));
    fields.push_back(std::to_string(matrix.columns));
    append_doubles(fields, matrix.values);
}

std::string canonical_parameter_contract(
    const epcsaft::NativeModel& model,
    const std::vector<EvaluatorParameter>& parameters
) {
    const auto coordinates = model.active_coordinates();
    if (coordinates.size() != parameters.size()) {
        throw std::invalid_argument("active parameter metadata does not match model");
    }
    const auto units = canonical_parameter_units(coordinates);
    const auto bases = canonical_parameter_bases(coordinates);
    std::string payload = "epcsaft.eos-owner-parameter-contract.v1|";
    append_field(payload, model.parameter_fingerprint());
    for (std::size_t index = 0; index < coordinates.size(); ++index) {
        const auto& parameter = parameters[index];
        if (parameter.identity != coordinates[index].identity
            || parameter.unit != units[index] || parameter.basis != bases[index]) {
            throw std::invalid_argument("active parameter metadata is not canonical");
        }
        append_field(payload, coordinates[index].identity);
        std::ostringstream value;
        value << std::scientific << std::setprecision(
            std::numeric_limits<double>::max_digits10
        ) << coordinates[index].value;
        append_field(payload, value.str());
        append_field(payload, units[index]);
        append_field(payload, bases[index]);
    }
    return "sha256:" + sha256_digest(payload);
}

void require_digest(const std::string& value, const char* name) {
    if (!sha256_fingerprint(value)) {
        throw std::invalid_argument(std::string(name) + " must be a sha256 fingerprint");
    }
}

std::string canonical_evaluator_fingerprint(
    const epcsaft::NativeModel& model,
    const GrepeRequest& request,
    const EquilibriumEvaluatorSidecar& sidecar
) {
    std::vector<std::string> fields{
        "equilibrium-owner-v1",
        model.parameter_fingerprint(),
        request.identity,
        sidecar.owner_kind,
        sidecar.source_spec_digest,
        sidecar.anchor_identity,
        sidecar.anchor_receipt_identity,
        sidecar.temperature_role,
        sidecar.pressure_role,
        sidecar.derivative_provenance,
        sidecar.topology_policy,
        sidecar.branch_policy,
        sidecar.conditioning_policy,
    };
    for (const auto& parameter : sidecar.active_parameters) {
        fields.insert(fields.end(), {
            parameter.identity, parameter.unit, parameter.basis
        });
    }
    for (const auto& phase : sidecar.phases) {
        fields.insert(fields.end(), {
            phase.identity, phase.role, phase.support,
            phase.finite ? "finite" : "incipient",
            phase.ideal_gas ? "ideal-gas" : "eos"
        });
        append_indices(fields, phase.support_indices);
    }
    fields.insert(fields.end(), model.component_ids().begin(), model.component_ids().end());
    for (const int charge : model.charge_numbers()) {
        fields.push_back(std::to_string(charge));
    }
    fields.insert(fields.end(), sidecar.reaction_identities.begin(),
        sidecar.reaction_identities.end());
    fields.insert(fields.end(), sidecar.reaction_phase_identities.begin(),
        sidecar.reaction_phase_identities.end());
    fields.insert(fields.end(), request.reaction.species_ids.begin(),
        request.reaction.species_ids.end());
    for (const int charge : request.reaction.charges) {
        fields.push_back(std::to_string(charge));
    }
    fields.push_back(request.reaction.eos_fingerprint);
    append_doubles(fields, request.reaction.molar_masses_kg_per_mol);
    append_matrix(fields, request.reaction.balance_matrix);
    append_matrix(fields, request.reaction.reaction_matrix);
    append_doubles(fields, request.reaction.feed_amounts);
    append_doubles(fields, request.reaction.ln_k);
    append_doubles(fields, request.reaction.ln_k_parameter_derivatives);
    for (std::size_t index = 0; index < sidecar.output_rows.size(); ++index) {
        const auto& row = sidecar.output_rows[index];
        fields.insert(fields.end(), {
            row.identity, row.unit, row.basis
        });
        const auto& output = request.outputs[index];
        fields.insert(fields.end(), {output.selector, std::to_string(output.phase_index)});
        append_doubles(fields, output.coefficients);
        append_doubles(fields, output.solvent_mass_coefficients_kg_per_mol);
    }
    fields.push_back(canonical_double(request.temperature.value));
    fields.push_back(canonical_double(request.pressure.value));
    append_doubles(fields, {
        request.temperature.bounds[0], request.temperature.bounds[1],
        request.pressure.bounds[0], request.pressure.bounds[1]
    });
    append_doubles(fields, request.temperature.starts);
    append_doubles(fields, request.pressure.starts);
    fields.push_back(canonical_double(request.continuation_pressure_pa));
    fields.push_back(canonical_double(request.trace_floor));
    fields.push_back(request.declared_starts_only ? "declared-starts-only" : "generated-starts");
    for (const auto& anchor : request.continuation_phases) {
        fields.push_back(anchor.identity);
        append_doubles(fields, anchor.mole_fractions);
        fields.push_back(canonical_double(anchor.molar_volume_m3_per_mol));
    }
    for (const auto& phase : request.phases) {
        fields.insert(fields.end(), {
            phase.model_reference_id,
            phase.ideal_gas ? "ideal-gas" : "eos",
            phase.finite ? "finite" : "incipient",
        });
        append_doubles(fields, {
            phase.packing_fraction_interval[0], phase.packing_fraction_interval[1]
        });
        if (phase.start.has_value()) {
            append_doubles(fields, phase.start->amounts_mol);
            fields.push_back(canonical_double(phase.start->molar_volume_m3_per_mol));
        }
    }
    return canonical_hash(fields);
}

EquilibriumEvaluatorSidecar make_sidecar(
    const epcsaft::NativeModel& model,
    const GrepeRequest& request,
    EquilibriumEvaluatorSetup setup,
    bool bubble_vle
) {
    const auto active = identities(model.active_coordinates());
    if (setup.active_parameters.size() != active.size()) {
        throw std::invalid_argument("active parameter metadata does not match model");
    }
    const auto coordinates = model.active_coordinates();
    const auto canonical_units = canonical_parameter_units(coordinates);
    const auto canonical_bases = canonical_parameter_bases(coordinates);
    for (std::size_t index = 0; index < active.size(); ++index) {
        if (setup.active_parameters[index].identity != active[index]
            || setup.active_parameters[index].unit != canonical_units[index]
            || setup.active_parameters[index].basis != canonical_bases[index]) {
            throw std::invalid_argument("active parameter metadata is not canonical");
        }
    }
    if (setup.output_rows.size() != request.outputs.size()
        || setup.output_rows.empty()) {
        throw std::invalid_argument("output metadata does not match request");
    }
    for (const auto& row : setup.output_rows) {
        if (row.identity.empty() || row.unit.empty() || row.basis.empty()) {
            throw std::invalid_argument("output metadata is incomplete");
        }
    }
    if (setup.anchor_identity.empty() || setup.anchor_receipt_identity.empty()) {
        throw std::invalid_argument("equilibrium evaluator identity metadata is incomplete");
    }
    if (setup.parameter_baseline_digest != model.parameter_fingerprint()) {
        throw std::invalid_argument("parameter baseline digest does not match model");
    }
    const std::string expected_source_digest =
        canonical_source_specification_fingerprint(
            request, !model.active_coordinates().empty()
        );
    if (setup.source_spec_digest != expected_source_digest) {
        throw std::invalid_argument("source specification digest does not match request");
    }
    if (!expected_source_digest.empty()) {
        require_digest(setup.source_spec_digest, "source specification digest");
    }
    std::vector<std::string> output_identities;
    output_identities.reserve(setup.output_rows.size());
    for (const auto& row : setup.output_rows) {
        output_identities.push_back(row.identity);
    }
    if (!unique_nonempty(active) || !unique_nonempty(output_identities)) {
        throw std::invalid_argument("evaluator identities are invalid");
    }
    EquilibriumEvaluatorSidecar result;
    result.owner_kind = setup.owner_kind.empty()
        ? (bubble_vle ? kBubbleVleKind : kHomogeneousKind)
        : setup.owner_kind;
    result.parameter_baseline_digest = std::move(setup.parameter_baseline_digest);
    result.source_spec_digest = std::move(setup.source_spec_digest);
    result.anchor_identity = std::move(setup.anchor_identity);
    result.anchor_receipt_identity = std::move(setup.anchor_receipt_identity);
    result.active_parameters = std::move(setup.active_parameters);
    result.output_rows = std::move(setup.output_rows);
    result.derivative_provenance = setup.derivative_provenance.empty()
        ? "equilibrium-solved-state-exact-total-derivative"
        : std::move(setup.derivative_provenance);
    result.phases.reserve(request.phases.size());
    for (const auto& phase : request.phases) {
        result.phases.push_back({
            phase.identity,
            phase.role == PhaseRole::Liquid ? "liquid" : "vapor",
            phase.support_name,
            phase.support_indices,
            phase.finite,
            phase.ideal_gas,
        });
    }
    result.temperature_role = "fixed";
    result.pressure_role = bubble_vle ? "solved" : "fixed";
    for (const auto& record : request.reaction.equilibrium_constant_records) {
        result.reaction_identities.push_back(record.reaction_id);
    }
    result.reaction_phase_identities = request.reaction_phase_ids;
    result.topology_policy = setup.topology_policy.empty()
        ? "declared-topology-only"
        : std::move(setup.topology_policy);
    result.branch_policy = setup.branch_policy.empty()
        ? "immutable-declared-local-branch"
        : std::move(setup.branch_policy);
    result.conditioning_policy = setup.conditioning_policy.empty()
        ? "typed-refusal-on-singular-or-ill-conditioned-sensitivity"
        : std::move(setup.conditioning_policy);
    result.parameter_contract_fingerprint = canonical_parameter_contract(
        model, result.active_parameters
    );
    result.evaluator_fingerprint = canonical_evaluator_fingerprint(
        model, request, result
    );
    return result;
}

interfaces::EvaluationResult map_result(
    const GrepeResult& result,
    const EquilibriumEvaluatorSidecar& sidecar,
    bool exact_jacobian
) {
    const auto& equilibrium = result.equilibrium;
    if (!equilibrium.accepted) {
        return failure(failure_code(result), equilibrium.failure_reason.empty()
            ? "declared equilibrium owner rejected the trial"
            : equilibrium.failure_reason);
    }
    const std::size_t output_count = sidecar.output_rows.size();
    const std::size_t parameter_count = sidecar.active_parameters.size();
    const bool parameter_order_matches = equilibrium.parameter_order.size() == parameter_count
        && std::equal(
            equilibrium.parameter_order.begin(), equilibrium.parameter_order.end(),
            sidecar.active_parameters.begin(),
            [](const std::string& identity, const EvaluatorParameter& parameter) {
                return identity == parameter.identity;
            }
        );
    if (!parameter_order_matches
        || equilibrium.output_values.size() != output_count
        || !finite(equilibrium.output_values)) {
        return failure("owner_contract_failure", "equilibrium output contract changed");
    }
    if (exact_jacobian) {
        const std::size_t width = equilibrium.input_order.size();
        if (equilibrium.sensitivity_status != "available"
            || equilibrium.output_input_derivatives.size() != output_count * width
            || !finite(equilibrium.output_input_derivatives)) {
            return failure(
                equilibrium.sensitivity_status == "unavailable"
                    ? "singular_sensitivity" : "owner_contract_failure",
                equilibrium.sensitivity_failure_reason.empty()
                    ? "exact equilibrium total derivatives are unavailable"
                    : equilibrium.sensitivity_failure_reason
            );
        }
    }
    interfaces::Evaluation evaluation;
    evaluation.predictions = equilibrium.output_values;
    if (!exact_jacobian) {
        return evaluation;
    }
    evaluation.jacobian.assign(output_count * parameter_count, 0.0);
    const std::size_t width = equilibrium.input_order.size();
    if (equilibrium.output_input_derivatives.size() != output_count * width) {
        return failure("owner_contract_failure", "equilibrium input derivative layout changed");
    }
    for (std::size_t output = 0; output < output_count; ++output) {
        for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
            const auto found = std::find(
                equilibrium.input_order.begin(), equilibrium.input_order.end(),
                sidecar.active_parameters[parameter].identity
            );
            if (found == equilibrium.input_order.end()) {
                return failure("owner_contract_failure", "equilibrium parameter order changed");
            }
            evaluation.jacobian[output * parameter_count + parameter] =
                equilibrium.output_input_derivatives[
                    output * width + static_cast<std::size_t>(found - equilibrium.input_order.begin())
                ];
        }
    }
    if (!finite(evaluation.jacobian)) {
        return failure("owner_contract_failure", "equilibrium Jacobian is nonfinite");
    }
    return evaluation;
}

} // namespace

EquilibriumOwnerEvaluator::EquilibriumOwnerEvaluator(
    epcsaft::NativeModel model,
    epcsaft::NativeModel incipient_model,
    GrepeRequest request,
    EquilibriumEvaluatorSidecar sidecar,
    bool bubble_vle
)
    : model_(std::move(model)), incipient_model_(std::move(incipient_model)),
      request_(std::move(request)), sidecar_(std::move(sidecar)), bubble_vle_(bubble_vle) {}

std::size_t EquilibriumOwnerEvaluator::parameter_count() const noexcept {
    return sidecar_.active_parameters.size();
}

std::size_t EquilibriumOwnerEvaluator::output_count() const noexcept {
    return sidecar_.output_rows.size();
}

const std::string& EquilibriumOwnerEvaluator::parameter_contract_fingerprint()
    const noexcept {
    return sidecar_.parameter_contract_fingerprint;
}

const std::string& EquilibriumOwnerEvaluator::evaluator_fingerprint()
    const noexcept {
    return sidecar_.evaluator_fingerprint;
}

interfaces::EvaluationResult EquilibriumOwnerEvaluator::evaluate(
    const std::vector<double>& physical_parameters,
    bool exact_jacobian
) const noexcept {
    try {
        if (physical_parameters.size() != parameter_count()) {
            return failure("parameter_order_mismatch", "active parameter dimensions changed");
        }
        if (!finite(physical_parameters)) {
            return failure("invalid_parameter", "active parameters must be finite");
        }
        const epcsaft::NativeModel model = model_at(model_, physical_parameters);
        GrepeResult result;
        if (bubble_vle_) {
            const std::vector<double> incipient_values = values_for_model(
                incipient_model_, sidecar_.active_parameters, physical_parameters
            );
            const epcsaft::NativeModel incipient_model = model_at(
                incipient_model_, incipient_values
            );
            result = solve_declared_bubble_vle(
                model, incipient_model, request_, exact_jacobian
            );
        } else {
            result = solve_declared_homogeneous_reactive(
                model, request_, exact_jacobian
            );
        }
        return map_result(result, sidecar_, exact_jacobian);
    } catch (const std::exception& error) {
        return invalid_result(error);
    } catch (...) {
        return failure("owner_failure", "equilibrium owner evaluation failed");
    }
}

std::shared_ptr<const EquilibriumOwnerEvaluator> compile_equilibrium_owner(
    const epcsaft::NativeModel& model,
    const epcsaft::NativeModel& incipient_model,
    GrepeRequest request,
    EquilibriumEvaluatorSetup setup
) {
    const bool bubble_vle = admits_declared_bubble_vle(request);
    const bool homogeneous = admits_declared_homogeneous_reactive(request);
    if (!bubble_vle && !homogeneous) {
        throw std::invalid_argument("only admitted declared reactive topologies are supported");
    }
    const std::string expected_kind = bubble_vle ? kBubbleVleKind : kHomogeneousKind;
    if (!setup.owner_kind.empty() && setup.owner_kind != expected_kind) {
        throw std::invalid_argument("equilibrium evaluator owner kind does not match topology");
    }
    request.require_output_derivatives = true;
    const EquilibriumEvaluatorSidecar sidecar = make_sidecar(
        model, request, std::move(setup), bubble_vle
    );
    if (sidecar.owner_kind != expected_kind) {
        throw std::invalid_argument("equilibrium evaluator owner kind is invalid");
    }
    const auto active = model.active_coordinates();
    const auto incipient_active = incipient_model.active_coordinates();
    const auto primary_ids = identities(active);
    const auto incipient_ids = identities(incipient_active);
    if (bubble_vle
        && model.parameter_fingerprint() != incipient_model.parameter_fingerprint()) {
        throw std::invalid_argument("bubble evaluator model fingerprints do not match");
    }
    if (!unique_nonempty(primary_ids) || !unique_nonempty(incipient_ids)) {
        throw std::invalid_argument("equilibrium active parameter identities are invalid");
    }
    for (const std::string& identity : incipient_ids) {
        if (std::find(primary_ids.begin(), primary_ids.end(), identity) == primary_ids.end()) {
            throw std::invalid_argument("incipient active parameter order changed");
        }
    }
    auto evaluator = std::shared_ptr<const EquilibriumOwnerEvaluator>(
        new EquilibriumOwnerEvaluator(
            model,
            incipient_model,
            std::move(request),
            sidecar,
            bubble_vle
        )
    );
    return evaluator;
}

} // namespace epcsaft_engine::equilibrium
