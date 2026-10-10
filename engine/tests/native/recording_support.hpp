#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <initializer_list>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace recording_support {

// execution_policy is read only by meson for equilibrium_checks (Automatic_fast or Explicit).
struct Declaration final {
    const char* test_id;
    const char* arguments;
    const char* execution_policy = nullptr;
};

using TestFunction = void (*)();

struct Case final {
    Declaration declaration;
    TestFunction function;
    const char* alias = nullptr;
};

class Failure final : public std::runtime_error {
public:
    explicit Failure(const std::string& message) : std::runtime_error(message) {}
};

inline void check(bool condition, const char* file, int line, const char* expression) {
    if (condition) return;
    std::fflush(stdout);
    std::fprintf(stderr, "%s:%d: failed: %s\n", file, line, expression);
    throw Failure(std::string(file) + ":" + std::to_string(line)
                  + ": failed: " + expression);
}

inline std::string number(double value) {
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.17g", value);
    return buffer;
}

inline std::string field(const char* name, double value) {
    return std::string(name) + ":" + number(value);
}

inline std::string field(const char* name, std::size_t value) {
    return std::string(name) + ":" + std::to_string(value);
}

template <class Integer, std::enable_if_t<std::is_integral_v<Integer>
                                          && !std::is_same_v<Integer, std::size_t>, int> = 0>
inline std::string field(const char* name, Integer value) {
    return std::string(name) + ":" + std::to_string(value);
}

inline std::string field(const char* name, const char* value) {
    return std::string(name) + ":" + value;
}

inline std::string field(const char* name, const std::string& value) {
    return std::string(name) + ":" + value;
}

inline std::string field(const char* name, const std::vector<double>& values) {
    std::string result = std::string(name) + ":";
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index != 0) result += ',';
        result += number(values[index]);
    }
    return result;
}

inline std::string field(const char* name, const std::vector<int>& values) {
    std::string result = std::string(name) + ":";
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index != 0) result += ',';
        result += std::to_string(values[index]);
    }
    return result;
}

inline std::string fields(std::initializer_list<std::string> values) {
    std::string result;
    for (const auto& value : values) {
        if (!result.empty()) result += ';';
        result += value;
    }
    return result;
}

// The key/value line is kept for the small native runner and for readable
// logs.  JSON records carry the same information without relying on a token
// delimiter, so claims, LaTeX labels and serialized inputs remain lossless.
inline std::string json_escape(const std::string& value) {
    std::string result;
    result.reserve(value.size() + 8);
    for (const unsigned char character : value) {
        switch (character) {
        case '"': result += "\\\""; break;
        case '\\': result += "\\\\"; break;
        case '\b': result += "\\b"; break;
        case '\f': result += "\\f"; break;
        case '\n': result += "\\n"; break;
        case '\r': result += "\\r"; break;
        case '\t': result += "\\t"; break;
        default:
            if (character < 0x20) {
                char buffer[7];
                std::snprintf(buffer, sizeof(buffer), "\\u%04x", character);
                result += buffer;
            } else {
                result += static_cast<char>(character);
            }
        }
    }
    return result;
}

inline std::string json_string(const char* value) {
    return std::string("\"") + json_escape(value == nullptr ? "" : value) + "\"";
}

inline std::string json_string(const std::string& value) {
    return std::string("\"") + json_escape(value) + "\"";
}

inline std::string json_number_or_null(double value) {
    return std::isfinite(value) ? number(value) : "null";
}

inline std::string json_scalar(const std::string& value) {
    char* end = nullptr;
    const double parsed = std::strtod(value.c_str(), &end);
    if (end != value.c_str() && *end == '\0' && std::isfinite(parsed))
        return number(parsed);
    return json_string(value);
}

inline std::string json_fields_object(const std::string& serialized) {
    std::string payload = "{";
    bool first = true;
    std::size_t start = 0;
    while (start <= serialized.size()) {
        const std::size_t end = serialized.find(';', start);
        const std::string entry = serialized.substr(
            start, end == std::string::npos ? std::string::npos : end - start
        );
        const std::size_t separator = entry.find(':');
        if (separator == std::string::npos) {
            if (!entry.empty()) {
                if (!first) payload += ',';
                payload += "\"raw\":[" + json_string(entry) + "]";
                first = false;
            }
        } else {
            const std::string key = entry.substr(0, separator);
            const std::string value = entry.substr(separator + 1);
            if (!first) payload += ',';
            payload += json_string(key) + ":";
            if (value.find(',') == std::string::npos) {
                payload += json_scalar(value);
            } else {
                payload += '[';
                std::size_t value_start = 0;
                bool value_first = true;
                while (value_start <= value.size()) {
                    const std::size_t value_end = value.find(',', value_start);
                    if (!value_first) payload += ',';
                    payload += json_scalar(value.substr(
                        value_start,
                        value_end == std::string::npos
                            ? std::string::npos : value_end - value_start
                    ));
                    value_first = false;
                    if (value_end == std::string::npos) break;
                    value_start = value_end + 1;
                }
                payload += ']';
            }
            first = false;
        }
        if (end == std::string::npos) break;
        start = end + 1;
    }
    payload += '}';
    return payload;
}

inline const char*& active_test_id() {
    static const char* value = nullptr;
    return value;
}

inline void print_json_compare(
    const char* comparison_id, const char* quantity, const char* unit,
    double model, std::optional<double> reference, const std::string& inputs, const char* basis,
    const char* method, const char* claim, const char* criterion,
    double atol, double rtol, bool asserted_verdict, const char* judgment,
    const char* comparison_group_id, const char* point_id, const char* quantity_latex,
    const char* inputs_latex, const char* relationship_latex, const char* coordinates,
    const char* criterion_rule = nullptr
) {
    std::string payload = "{\"schema_version\":1";
    payload += ",\"comparison_id\":" + json_string(comparison_id);
    payload += ",\"quantity\":" + json_string(quantity);
    payload += ",\"unit\":" + json_string(unit);
    payload += ",\"model\":" + json_number_or_null(model);
    if (reference.has_value())
        payload += ",\"reference\":" + json_number_or_null(*reference);
    payload += ",\"inputs\":" + json_fields_object(inputs);
    payload += ",\"basis\":" + json_string(basis);
    payload += ",\"method\":" + json_string(method);
    payload += ",\"claim\":" + json_string(claim);
    payload += ",\"judgment\":" + json_string(judgment);
    const bool finite_operands = std::isfinite(model)
        && (!reference.has_value() || std::isfinite(*reference));
    payload += ",\"evidence_state\":"
        + json_string(finite_operands ? "Complete" : "Unavailable");
    if (!finite_operands)
        payload += ",\"unavailable_reason\":\"nonfinite comparison operand\"";
    if (criterion != nullptr) {
        payload += ",\"criterion\":" + json_string(criterion);
        payload += ",\"rule\":"
            + json_string(criterion_rule == nullptr ? criterion : criterion_rule);
        if (std::string(criterion) != "custom") {
            payload += ",\"atol\":" + json_number_or_null(atol);
            payload += ",\"rtol\":" + json_number_or_null(rtol);
        }
        payload += ",\"verdict\":";
        payload += asserted_verdict ? "true" : "false";
    }
    if (comparison_group_id != nullptr)
        payload += ",\"comparison_group_id\":" + json_string(comparison_group_id);
    if (point_id != nullptr)
        payload += ",\"point_id\":" + json_string(point_id);
    if (quantity_latex != nullptr)
        payload += ",\"quantity_latex\":" + json_string(quantity_latex);
    if (relationship_latex != nullptr)
        payload += ",\"relationship_latex\":" + json_string(relationship_latex);
    if (inputs_latex != nullptr)
        payload += ",\"inputs_latex\":" + json_string(inputs_latex);
    if (coordinates != nullptr)
        payload += ",\"coordinates\":" + json_fields_object(coordinates);
    payload += '}';
    std::printf("EPCSAFT_COMPARE_JSON %s\n", payload.c_str());
    std::fflush(stdout);
}

inline bool scalar_verdict(double model, double reference, const char* criterion,
                           double atol, double rtol) {
    const std::string rule = criterion == nullptr ? std::string() : std::string(criterion);
    const bool known = rule == "symmetric_close"
        || rule == "reference_scale"
        || rule == "strict_absolute"
        || rule == "absolute"
        || rule == "strict_symmetric_close"
        || rule == "strict_reference_scale"
        || rule == "strict_reference_or_one"
        || rule == "reference_or_one";
    if (!known) throw std::invalid_argument("unknown scalar criterion: " + rule);
    if (!std::isfinite(model) || !std::isfinite(reference)
        || !std::isfinite(atol) || !std::isfinite(rtol)
        || atol < 0.0 || rtol < 0.0) {
        return false;
    }
    const double error = std::abs(model - reference);
    if (rule == "symmetric_close") {
        return error <= std::max(atol, rtol * std::max(std::abs(model), std::abs(reference)));
    }
    if (rule == "reference_scale") {
        return error <= atol + rtol * std::abs(reference);
    }
    if (rule == "strict_absolute") {
        return error < atol;
    }
    if (rule == "absolute") {
        return error <= atol;
    }
    if (rule == "strict_symmetric_close") {
        return error < std::max(atol, rtol * std::max(std::abs(model), std::abs(reference)));
    }
    if (rule == "strict_reference_scale") {
        return error < atol + rtol * std::abs(reference);
    }
    if (rule == "strict_reference_or_one") {
        return error < atol + rtol * std::max(1.0, std::abs(reference));
    }
    return error <= atol + rtol * std::max(1.0, std::abs(reference));
}

inline void print(const Declaration& declaration, const char* alias) {
    std::printf("EPCSAFT_DECL schema_version=1 test_id=%s arguments=%s", declaration.test_id, declaration.arguments);
    if (declaration.execution_policy != nullptr)
        std::printf(" execution_policy=%s", declaration.execution_policy);
    if (alias != nullptr) std::printf(" alias=%s", alias);
    std::putchar('\n');
    std::fflush(stdout);
}

inline bool compare(
    const char* comparison_id, const char* quantity, const char* unit,
    double model, double reference, const std::string& inputs, const char* basis,
    const char* method, const char* claim, const char* criterion,
    double atol, double rtol, bool asserted_verdict, const char* comparison_group_id = nullptr,
    const char* point_id = nullptr, const char* quantity_latex = nullptr,
    const char* inputs_latex = nullptr, const char* relationship_latex = nullptr,
    const char* coordinates = nullptr
) {
    const bool verdict = scalar_verdict(model, reference, criterion, atol, rtol);
    print_json_compare(
        comparison_id, quantity, unit, model, reference, inputs, basis, method, claim,
        criterion, atol, rtol, verdict, verdict ? "Criterion satisfied" : "Criterion failed",
        comparison_group_id, point_id, quantity_latex, inputs_latex, relationship_latex,
        coordinates
    );
    if (verdict != asserted_verdict)
        throw std::logic_error(
            std::string("recorded criterion disagrees with assertion: ") + comparison_id);
    if (!verdict) throw Failure(std::string("criterion failed: ") + comparison_id);
    return verdict;
}

inline bool compare_custom(
    const char* comparison_id, const char* quantity, const char* unit,
    double model, double reference, const std::string& inputs, const char* basis,
    const char* method, const char* claim, const char* rule, bool asserted_verdict,
    const char* comparison_group_id = nullptr, const char* point_id = nullptr,
    const char* quantity_latex = nullptr, const char* inputs_latex = nullptr,
    const char* relationship_latex = nullptr, const char* coordinates = nullptr
) {
    print_json_compare(
        comparison_id, quantity, unit, model, reference, inputs, basis, method, claim,
        "custom", 0.0, 0.0, asserted_verdict,
        asserted_verdict ? "Criterion satisfied" : "Criterion failed",
        comparison_group_id, point_id, quantity_latex, inputs_latex, relationship_latex,
        coordinates, rule
    );
    if (!asserted_verdict) throw Failure(std::string("criterion failed: ") + comparison_id);
    return asserted_verdict;
}

inline void observe(
    const char* comparison_id, const char* quantity, const char* unit,
    double model, std::optional<double> reference, const std::string& inputs, const char* basis,
    const char* method, const char* claim, const char* comparison_group_id = nullptr,
    const char* point_id = nullptr, const char* quantity_latex = nullptr,
    const char* inputs_latex = nullptr, const char* relationship_latex = nullptr,
    const char* coordinates = nullptr
) {
    print_json_compare(
        comparison_id, quantity, unit, model, reference, inputs, basis, method, claim,
        nullptr, 0.0, 0.0, false, "Observation only", comparison_group_id, point_id,
        quantity_latex, inputs_latex, relationship_latex, coordinates
    );
}

inline bool listing(int argc, char** argv) {
    return argc == 2 && argv != nullptr && argv[1] != nullptr
        && std::string(argv[1]) == "--list";
}

inline void case_event(const char* test_id, const char* execution, const char* reason = nullptr) {
    std::printf("EPCSAFT_CASE schema_version=1 test_id=%s execution=%s", test_id, execution);
    if (reason != nullptr) std::printf(" outcome_reason=%s", reason);
    std::putchar('\n');
    std::fflush(stdout);
}

inline int run_case(const Case& item) {
    case_event(item.declaration.test_id, "Running");
    active_test_id() = item.declaration.test_id;
    try {
        item.function();
        case_event(item.declaration.test_id, "Passed");
        active_test_id() = nullptr;
        return 0;
    } catch (const Failure& error) {
        case_event(item.declaration.test_id, "Failed", error.what());
        active_test_id() = nullptr;
        return 1;
    } catch (const std::exception& error) {
        case_event(item.declaration.test_id, "Error", error.what());
        active_test_id() = nullptr;
        return 1;
    } catch (...) {
        case_event(item.declaration.test_id, "Error", "unknown native exception");
        active_test_id() = nullptr;
        return 1;
    }
}

inline bool selected(const Case& item, const std::string& argument) {
    return argument == item.declaration.test_id
        || argument == (item.declaration.arguments == nullptr ? "" : item.declaration.arguments)
        || (item.alias != nullptr && argument == item.alias);
}

inline int run(int argc, char** argv, std::initializer_list<Case> cases) {
    if (listing(argc, argv)) {
        for (const auto& item : cases) print(item.declaration, item.alias);
        return 0;
    }
    if (argc > 2) {
        std::fprintf(stderr, "expected one case selector; --list shows test ids, arguments and aliases\n");
        return 2;
    }
    int status = 0;
    if (argc == 2) {
        // A test id or argument selects one case; an alias shared by several cases runs them all.
        const std::string argument = argv[1] == nullptr ? "" : argv[1];
        bool any = false;
        for (const auto& item : cases) {
            if (!selected(item, argument)) continue;
            any = true;
            status |= run_case(item);
        }
        if (!any) std::fprintf(stderr, "no case matches %s; --list shows test ids, arguments and aliases\n", argument.c_str());
        return any ? status : 2;
    }
    for (const auto& item : cases) status |= run_case(item);
    return status;
}

}  // namespace recording_support
