#if defined(_WIN32)
#define EPC_DEFINE_POLICY
#endif
#include "core/trace.hpp"
#include "core/issue235_prototype.hpp"

namespace epcsaft::issue235 {
thread_local Row row;
extern "C" EPC_API void ep235_prototype_begin(unsigned mode) {
    const auto generation = row.generation + 1;
    row = {};
    row.generation = generation;
    row.mode = mode;
}
extern "C" EPC_API const std::uint64_t* ep235_prototype_end() {
    row.mode = 0;
    return row.counts.data();
}
}

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
