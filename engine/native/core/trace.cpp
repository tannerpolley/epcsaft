#include "core/trace.hpp"

namespace epcsaft::core {
thread_local Trace* active_trace = nullptr;
Trace*& current_trace() noexcept { return active_trace; }
}  // namespace epcsaft::core
