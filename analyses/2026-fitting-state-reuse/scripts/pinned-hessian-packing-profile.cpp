// Instrument the immutable wheel; every calculation still calls its original owner.
#include "core/trace.hpp"
#include "equilibrium/problem.hpp"
#include "properties/helmholtz_model.hpp"
#include <cstdlib>
#include <dlfcn.h>
#include <stdexcept>

namespace {
void* original(const char* symbol) {
    static void* module = dlopen(std::getenv("EPCSAFT_PROFILE_CORE"), RTLD_NOW | RTLD_NOLOAD);
    void* function = module ? dlsym(module, symbol) : nullptr;
    if (!function) throw std::runtime_error("immutable wheel profiling symbol unavailable");
    return function;
}
// The wheel hides its active_trace TLS; use its constructor and the unchanged inline destructor.
class NativeScope final {
public:
    explicit NativeScope(const char* name) {
        using Constructor = void (*)(epcsaft::core::Scope*, const char*);
        static auto construct = reinterpret_cast<Constructor>(original("_ZN7epcsaft4core5ScopeC1EPKc"));
        construct(reinterpret_cast<epcsaft::core::Scope*>(storage_), name);
    }
    ~NativeScope() { reinterpret_cast<epcsaft::core::Scope*>(storage_)->~Scope(); }
private:
    alignas(epcsaft::core::Scope) unsigned char storage_[sizeof(epcsaft::core::Scope)];
};
}

namespace epcsaft {
double HelmholtzModel::packing_fraction_coefficient(double T, const std::vector<double>& x) const {
    using Function = double (*)(const HelmholtzModel*, double, const std::vector<double>&);
    static auto function = reinterpret_cast<Function>(original(
        "_ZNK7epcsaft14HelmholtzModel28packing_fraction_coefficientEdRKSt6vectorIdSaIdEE"));
    NativeScope scope("packing_fraction_coefficient");
    return function(this, T, x);
}
}

namespace epcsaft::equilibrium {
std::vector<double> CompiledSystem::constraint_hessian_at(
    const PhasePoint& point, const std::vector<double>& multipliers) const {
    using Function = std::vector<double> (*)(const CompiledSystem*, const PhasePoint&, const std::vector<double>&);
    static auto function = reinterpret_cast<Function>(original(
        "_ZNK7epcsaft11equilibrium14CompiledSystem21constraint_hessian_atERKNS0_10PhasePointERKSt6vectorIdSaIdEE"));
    NativeScope scope(pinned_count ? "pinned_hessian" : "constraint_hessian");
    return function(this, point, multipliers);
}
}
