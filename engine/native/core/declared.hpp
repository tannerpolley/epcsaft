#pragma once

#include <optional>
#include <stdexcept>
#include <string>

namespace epcsaft::core {

// A missing or physically inadmissible declaration.
struct __attribute__((visibility("default"))) DeclarationError final : std::invalid_argument {
    using std::invalid_argument::invalid_argument;
};

template <typename T>
class Declared final {
public:
    constexpr Declared() = default;
    constexpr Declared(T value) : value_(value) {}

    constexpr bool declared() const noexcept { return value_.has_value(); }
    constexpr bool is(T value) const noexcept {
        return value_.has_value() && *value_ == value;
    }
    const T& require(const char* name) const {
        if (!value_.has_value()) throw DeclarationError(std::string("missing declaration: ") + name);
        return *value_;
    }

private:
    std::optional<T> value_;
};

}  // namespace epcsaft::core
