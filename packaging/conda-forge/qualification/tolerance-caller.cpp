#include "core/numerical_policy.hpp"
#include "core/declared.hpp"
#include "equilibrium/problem.hpp"
#include "core/trace.hpp"
#include <iostream>
#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif
int main() {
#ifdef _WIN32
    for (const wchar_t* name : {L"epcsaft_native.dll", L"epcsaft_ipopt_3119.dll", L"dmumps.dll", L"openblas.dll"}) {
        char path[32768];
        const auto module = GetModuleHandleW(name);
        if (!module || !GetModuleFileNameA(module, path, 32768)) return 1;
        std::cout << "loaded_module=" << path << "\n";
    }
#else
    // The image defining one exported function each: the native library, the bundled Ipopt, MUMPS and BLAS.
    // A function, not data: an executable may hold a copy of the library's data.
    const void* symbols[] = {reinterpret_cast<const void*>(&epcsaft::core::current_trace), dlsym(RTLD_DEFAULT, "CreateIpoptProblem"),
                             dlsym(RTLD_DEFAULT, "dmumps_c"), dlsym(RTLD_DEFAULT, "dgemm_")};
    for (const void* symbol : symbols) {
        Dl_info image{};
        if (!symbol || !dladdr(symbol, &image) || !image.dli_fname) return 1;
        std::cout << "loaded_module=" << image.dli_fname << "\n";
    }
#endif
    bool caught=false;
    using namespace epcsaft::equilibrium;IdealMoleFraction model({"a"},{0});EquilibriumProblem request;request.temperature=FreeValue::fixed(300);request.pressure=FreeValue::fixed(1e5);request.feed.values={1};PhaseSpec phase;phase.name="L";phase.kind=PhaseKind::Liquid;phase.support_declared=true;request.phases={phase};try{compile_problem(model,request);}catch(const epcsaft::core::DeclarationError& error){std::cout<<"DeclarationError caught across native library: "<<error.what()<<'\n';caught=true;}catch(const std::exception& error){std::cerr<<"unexpected type: "<<error.what()<<'\n';return 1;}
    if(!caught)return 1;
    const auto* imported = &epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL;
    const epcsaft::policy::Tolerance* table = nullptr;
    for (const auto* row : epcsaft::policy::ENTRIES) if (row->name == imported->name) table = row;
    std::cout << "shared tolerance imported=" << imported << " table=" << table << "\n";
    return imported == table ? 0 : 1;
}
