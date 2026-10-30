// HS071 and its analytic derivatives are supplied by the locked Ipopt archive.
#include "hs071_nlp.hpp"
#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include "IpIpoptData.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpTNLPAdapter.hpp"
#include "IpCompoundVector.hpp"
#include "IpIteratesVector.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <dlfcn.h>
#endif

class InternalHS071 : public HS071_NLP {
public:
  bool observed = false;
  bool intermediate_callback(AlgorithmMode, Index iter, Number, Number, Number,
      Number, Number, Number, Number, Number, Index,
      const IpoptData* data, IpoptCalculatedQuantities* cq) override {
    IpoptNLP& active = *cq->GetIpoptNLP();
    auto* original = dynamic_cast<OrigIpoptNLP*>(&active);
    auto* adapter = original ? dynamic_cast<TNLPAdapter*>(GetRawPtr(original->nlp())) : nullptr;
    if (!adapter || active.NLP_scaling()->have_x_scaling()) return false;
    const Vector* current = GetRawPtr(data->curr()->x());
    if (const auto* compound = dynamic_cast<const CompoundVector*>(current)) current = GetRawPtr(compound->GetComp(0));
    SmartPtr<const Vector> lower = active.x_L(), upper = active.x_U();
    double x[4], lo[4], hi[4];
    adapter->ResortX(*current, x);
    adapter->ResortBnds(*lower, lo, *upper, hi);
    for (int i = 0; i < 4; ++i) {
      if (!std::isfinite(x[i]) || !std::isfinite(lo[i]) || !std::isfinite(hi[i])) return false;
      std::cout << "internal," << iter << ',' << i << ',' << x[i] << ',' << lo[i] << ',' << hi[i] << '\n';
    }
    observed = true;
    return true;
  }
  void finalize_solution(SolverReturn status, Index n, const Number* x,
      const Number* zl, const Number* zu, Index m, const Number* g,
      const Number* lambda, Number obj, const IpoptData* data,
      IpoptCalculatedQuantities* cq) override {
    HS071_NLP::finalize_solution(status, n, x, zl, zu, m, g, lambda, obj, data, cq);
    for (int i = 0; i < m; ++i) std::cout << "lambda[" << i << "] = " << lambda[i] << '\n';
  }
};
int main() {
  static_assert(IPOPT_VERSION_MAJOR == 3 && IPOPT_VERSION_MINOR == 11 && IPOPT_VERSION_RELEASE == 9, "locked Ipopt required");
  std::cout << std::setprecision(17) << "version=" << IPOPT_VERSION << "\nlinear_solver=mumps\n";
#ifdef _WIN32
  HMODULE modules[1024]; DWORD needed;
  if (!EnumProcessModules(GetCurrentProcess(), modules, sizeof(modules), &needed) || needed > sizeof(modules)) return 3;
  for (unsigned i = 0; i < needed / sizeof(HMODULE); ++i) {
    char path[MAX_PATH];
    if (!GetModuleFileNameA(modules[i], path, MAX_PATH)) return 3;
    std::cout << "module=" << path << '\n';
  }
#else
  Dl_info info;
  if (!dladdr(reinterpret_cast<void*>(&IpoptApplicationFactory), &info)) return 3;
  std::cout << "module=" << info.dli_fname << '\n';
#endif
  SmartPtr<InternalHS071> nlp = new InternalHS071;
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
  app->Options()->SetNumericValue("tol", 1e-7);
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetStringValue("linear_solver", "mumps");
  app->Options()->SetStringValue("option_file_name", "");
  if (app->Initialize() != Solve_Succeeded) return 2;
  Number relaxation;
  app->Options()->GetNumericValue("bound_relax_factor", relaxation, "");
  std::cout << "bound_relax_factor=" << relaxation << '\n';
  auto status = app->OptimizeTNLP(nlp);
  auto stats = app->Statistics();
  Number dual, constr, complementarity, kkt;
  stats->Infeasibilities(dual, constr, complementarity, kkt);
  std::cout << "diagnostics=" << dual << ',' << constr << ',' << complementarity << ',' << kkt << '\n';
  stats->ScaledInfeasibilities(dual, constr, complementarity, kkt);
  std::cout << "scaled_diagnostics=" << dual << ',' << constr << ',' << complementarity << ',' << kkt << '\n';
  Index f, g, df, dg, h;
  stats->NumberOfEvaluations(f, g, df, dg, h);
  std::cout << "work=" << stats->IterationCount() << ',' << f << ',' << g << ',' << df << ',' << dg << ',' << h << '\n';
  std::cout << "status=" << int(status) << '\n';
  return status == Solve_Succeeded && nlp->observed && std::isfinite(kkt) && kkt <= 1e-7 ? 0 : 1;
}
