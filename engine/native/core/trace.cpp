#if defined(_WIN32)
#define EPC_DEFINE_POLICY
#endif
#include "core/trace.hpp"

namespace epcsaft::core {
thread_local Trace* active_trace = nullptr;
Trace*& current_trace() noexcept { return active_trace; }
CalculationWork*& CalculationWork::current() noexcept {
    static thread_local CalculationWork* active = nullptr;
    return active;
}
AssociationObservations*& AssociationObservations::current() noexcept {
    static thread_local AssociationObservations* active = nullptr;
    return active;
}

}  // namespace epcsaft::core
