#pragma once

#include <optional>
#include <stdexcept>
#include <string>

namespace epcsaft::core {

struct MissingDeclaration final : std::invalid_argument {
    explicit MissingDeclaration(const char* name)
        : std::invalid_argument(std::string("missing declaration: ") + name) {}
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
        if (!value_.has_value()) throw MissingDeclaration(name);
        return *value_;
    }

private:
    std::optional<T> value_;
};

}  // namespace epcsaft::core
