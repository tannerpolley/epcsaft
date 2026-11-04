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
#include <vector>

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

struct ConditioningObservation {
    std::string operation{}, basis{}, norm{}, method{}, cause = "not_evaluated", point_kind = "last_linearization", criterion{};
    std::size_t dimension = 0, rank = 0, precision_bits = std::numeric_limits<double>::digits, iteration = 0;
    std::optional<std::size_t> attempt{}, stage{}, phase{}, component{};
    double threshold = 0, condition = INFINITY, matrix_norm = NAN, error = NAN, criterion_value = NAN, export_displacement = NAN;
    std::vector<double> coordinates{}, row_scale{}, column_scale{}, controls{};
    bool available = false;
};
struct NumericalRank final {
    std::string basis;
    std::size_t dimension = 0, rank = 0;
    double threshold = 0;
    ConditioningObservation observation{};
    std::shared_ptr<ConditioningObservation> first_conditioning_warning{};
};
inline void observe_conditioning(NumericalRank& record, ConditioningObservation incoming,
    const std::vector<double>& point = {}) {
    if (incoming.error == INFINITY) incoming.cause = "nonfinite_projection";
    incoming.basis = record.basis; incoming.rank = record.rank; incoming.dimension = record.dimension;
    incoming.threshold = record.threshold;
    if (!record.first_conditioning_warning && (incoming.cause == "rank_deficient"
        || incoming.cause == "nonfinite_inverse" || incoming.cause == "association_pivot_refusal"
        || (incoming.error == INFINITY || incoming.error > incoming.criterion_value))) {
        record.first_conditioning_warning = std::make_shared<ConditioningObservation>(incoming);
        record.first_conditioning_warning->coordinates = point;
    }
    record.observation = std::move(incoming);
}

struct Convergence final {
    double residual = 0.0;
    const policy::Tolerance* criterion = nullptr;
    std::optional<double> floor, amplification, contraction;
    std::size_t iterations = 0;
    std::size_t evaluations = 0;
    Stop stop = Stop::stalled;
    // Owns a caller-valued criterion; static policy producers leave it empty.
    std::shared_ptr<const policy::Tolerance> criterion_owner;
    std::string cause;
    std::shared_ptr<NumericalRank> rank;
    std::optional<double> backward_error;
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
    if (!record.cause.empty()) result += "; " + record.cause;
    if (record.rank) result += "; " + record.rank->basis + " rank " + std::to_string(record.rank->rank)
        + "/" + std::to_string(record.rank->dimension) + "; threshold " + optional_number(record.rank->threshold);
    if (record.backward_error) result += "; backward error " + optional_number(record.backward_error);
    return result;
}

// ALG-ACTIONS verification only: an exact derivative A against D = (f+ - f-)/(2h) of its value (#976).
struct DifferenceCheck final {
    std::string label, status = "nonfinite", cause;  // cause: why a neighbor was unavailable
    double exact = NAN, difference = NAN, step = NAN, roundoff = NAN, ratio = NAN;
};
inline void admit_difference_parameters(double step, double tau, double kappa) {
    if (!(std::isfinite(step) && step > 0.0 && std::isfinite(tau) && tau > 0.0 && std::isfinite(kappa) && kappa > 0.0))
        throw std::invalid_argument("difference step, tau and kappa must be finite and positive");
}
// Roundoff r = (e+ + e-)/(2h) from each value's evaluation floor e. Passes when |D - A| <= tau |A| + kappa r; a zero
// allowance passes only a zero error (ratio 0, else +inf); any nonfinite operand is "nonfinite" with no ratio.
inline DifferenceCheck difference_check(std::string label, double exact, double plus, double minus, double floor_plus,
                                        double floor_minus, double step, double tau, double kappa) {
    // A finite difference implies finite plus and minus.
    const double d = (plus - minus) / (2.0 * step), r = (floor_plus + floor_minus) / (2.0 * step);
    DifferenceCheck c{std::move(label), "nonfinite", {}, exact, d, step, r};
    if (!(std::isfinite(exact) && std::isfinite(d) && std::isfinite(r))) return c;
    const double error = std::abs(c.difference - exact), allowance = tau * std::abs(exact) + kappa * c.roundoff;
    c.ratio = allowance > 0.0 ? error / allowance : error == 0.0 ? 0.0 : INFINITY;
    c.status = c.ratio <= 1.0 ? "passed" : "failed";
    return c;
}

class EPC_EXCEPTION ConvergenceFailure : public std::runtime_error {
public:
    explicit ConvergenceFailure(const Convergence& record)
        : std::runtime_error(describe(record)), record(record) {}
    ConvergenceFailure(const Convergence& record, const std::optional<Convergence>& first)
        : ConvergenceFailure(record) { first_numerical_failure = first; }

    Convergence record;
    std::optional<Convergence> first_numerical_failure;
};

}  // namespace epcsaft::core
