#pragma once

#include "core/numerical_policy.hpp"

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <deque>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <cstdint>

namespace epcsaft::core {

enum class Stop : std::uint8_t {
    converged,
    converged_to_floor,
    stalled,
    budget_exhausted,
    left_domain,
    singular,
    infeasible,
    nonfinite,
};

struct Convergence final {
    double residual = 0.0;
    const policy::Tolerance* criterion = nullptr;
    std::optional<double> floor;
    std::optional<double> amplification;
    std::optional<double> contraction;
    std::size_t iterations = 0;
    std::size_t evaluations = 0;
    Stop stop = Stop::stalled;
    // Owns a caller-valued criterion; static policy producers leave it empty.
    std::shared_ptr<const policy::Tolerance> criterion_owner;
};

inline void stopping_criterion(Convergence& record, const policy::SolverOptions& options, std::size_t i) {
    record.criterion_owner = options.owners.at(i);
    record.criterion = options.criterion(i);
}

inline Stop classify(
    double residual,
    const policy::Tolerance& criterion,
    std::optional<double> floor,
    std::optional<double> contraction,
    bool budget_hit,
    std::optional<double> amplification = std::nullopt
) noexcept {
    if (!std::isfinite(residual) ||
        (floor.has_value() && !std::isfinite(*floor)) ||
        (contraction.has_value() && !std::isfinite(*contraction)) ||
        (amplification.has_value() && !std::isfinite(*amplification))) {
        return Stop::nonfinite;
    }
    if (residual <= criterion.value) return Stop::converged;
    if (floor.has_value() && residual <= policy::FLOOR_MARGIN.value * *floor) {
        return Stop::converged_to_floor;
    }
    if (budget_hit && (!contraction.has_value() || *contraction < 1.0)) {
        return Stop::budget_exhausted;
    }
    return Stop::stalled;
}

class ResidualHistory final {
public:
    void push(double residual) {
        if (!std::isfinite(residual)) throw std::invalid_argument("nonfinite residual history");
        values_.push_back(residual);
        const auto limit = static_cast<std::size_t>(policy::STALL_WINDOW.value) + 1;
        if (values_.size() > limit) values_.pop_front();
    }

    std::optional<double> contraction() const {
        const auto window = static_cast<std::size_t>(policy::STALL_WINDOW.value);
        if (values_.size() <= window) return std::nullopt;
        const double first = values_.front();
        const double last = values_.back();
        if (!(first > 0.0) || !std::isfinite(first) || !std::isfinite(last)) {
            return std::nullopt;
        }
        return std::pow(last / first, 1.0 / static_cast<double>(window));
    }

private:
    std::deque<double> values_;
};

inline const char* stop_name(Stop stop) noexcept {
    switch (stop) {
    case Stop::converged: return "converged";
    case Stop::converged_to_floor: return "converged_to_floor";
    case Stop::stalled: return "stalled";
    case Stop::budget_exhausted: return "budget_exhausted";
    case Stop::left_domain: return "left_domain";
    case Stop::singular: return "singular";
    case Stop::infeasible: return "infeasible";
    case Stop::nonfinite: return "nonfinite";
    }
    return "unknown";
}

inline std::string optional_number(const std::optional<double>& value) {
    if (!value.has_value()) return "absent";
    char buffer[64] = {};
    std::snprintf(buffer, sizeof(buffer), "%.17g", *value);
    return buffer;
}

inline std::string describe(const Convergence& record) {
    const char* criterion_name = record.criterion == nullptr
        ? "absent" : record.criterion->name;
    const double criterion_value = record.criterion == nullptr
        ? 0.0 : record.criterion->value;
    const std::string contraction = optional_number(record.contraction);
    const std::string floor = optional_number(record.floor);
    const std::string amplification = optional_number(record.amplification);
    std::string result(512, '\0');
    const int size = std::snprintf(
        result.data(), result.size(),
        "%s: residual %.17g against %s %.17g after %zu iterations, contraction %s, floor %s, amplification %s",
        stop_name(record.stop), record.residual, criterion_name, criterion_value,
        record.iterations, contraction.c_str(), floor.c_str(), amplification.c_str()
    );
    if (size < 0) return {};
    result.resize(static_cast<std::size_t>(size));
    return result;
}

class __attribute__((visibility("default"))) ConvergenceFailure : public std::runtime_error {
public:
    explicit ConvergenceFailure(const Convergence& record)
        : std::runtime_error(describe(record)), record(record) {}

    Convergence record;
};

}  // namespace epcsaft::core
