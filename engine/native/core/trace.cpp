#include "core/trace.hpp"

namespace epcsaft::core {
thread_local Trace* active_trace = nullptr;
Trace*& current_trace() noexcept { return active_trace; }

ConvergenceFailure::ConvergenceFailure(const Convergence& value, std::optional<Convergence>& first) : ConvergenceFailure(value) {
    if (!observe_diagnostic([&] { first_numerical_failure = first; })) {
        first->incomplete = true;
        first_numerical_failure = std::move(first);
    }
}

DeclarationError::DeclarationError(const std::string& message, Convergence& value) : std::invalid_argument(message) {
    if (!observe_diagnostic([&] { record = value; })) {
        value.incomplete = true;
        record = std::move(value);
    }
}
}  // namespace epcsaft::core
