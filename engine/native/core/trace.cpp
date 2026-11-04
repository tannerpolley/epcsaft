#if defined(_WIN32)
#define EPC_DEFINE_POLICY
#endif
#include "core/trace.hpp"

namespace epcsaft::core {
thread_local Trace* active_trace = nullptr;
Trace*& current_trace() noexcept { return active_trace; }
AssociationObservations*& AssociationObservations::current() noexcept {
    static thread_local AssociationObservations* active = nullptr;
    return active;
}

}  // namespace epcsaft::core
