#pragma once
#include "core/convergence.hpp"
#include "equilibrium/problem.hpp"
#include "core/trace.hpp"
#if __has_include(<HSLLoader.h>)
#include <HSLLoader.h>
#endif
#include <IpIpoptApplication.hpp>
#include <IpTNLP.hpp>
#include <algorithm>
#include <cmath>
#include <map>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include "core/numerical_policy.hpp"
#include <stdexcept>
#include <string>
#include <utility>

namespace epcsaft::equilibrium {
inline const char* solver_status_name(Ipopt::SolverReturn status) noexcept {
    switch (status) {
    case Ipopt::SUCCESS: return "SUCCESS";
    case Ipopt::MAXITER_EXCEEDED: return "MAXITER_EXCEEDED";
    case Ipopt::CPUTIME_EXCEEDED: return "CPUTIME_EXCEEDED";
    case Ipopt::STOP_AT_TINY_STEP: return "STOP_AT_TINY_STEP";
    case Ipopt::STOP_AT_ACCEPTABLE_POINT: return "STOP_AT_ACCEPTABLE_POINT";
    case Ipopt::LOCAL_INFEASIBILITY: return "LOCAL_INFEASIBILITY";
    case Ipopt::USER_REQUESTED_STOP: return "USER_REQUESTED_STOP";
    case Ipopt::FEASIBLE_POINT_FOUND: return "FEASIBLE_POINT_FOUND";
    case Ipopt::DIVERGING_ITERATES: return "DIVERGING_ITERATES";
    case Ipopt::RESTORATION_FAILURE: return "RESTORATION_FAILURE";
    case Ipopt::ERROR_IN_STEP_COMPUTATION: return "ERROR_IN_STEP_COMPUTATION";
    case Ipopt::INVALID_NUMBER_DETECTED: return "INVALID_NUMBER_DETECTED";
    case Ipopt::TOO_FEW_DEGREES_OF_FREEDOM: return "TOO_FEW_DEGREES_OF_FREEDOM";
    case Ipopt::INVALID_OPTION: return "INVALID_OPTION";
    case Ipopt::OUT_OF_MEMORY: return "OUT_OF_MEMORY";
    case Ipopt::INTERNAL_ERROR: return "INTERNAL_ERROR";
    case Ipopt::UNASSIGNED: return "UNASSIGNED";
    }
    return "UNKNOWN";
}

// The fatal Ipopt application returns (IpReturnCodes_inc.h) by name; any other code by value.
inline std::string application_return_name(Ipopt::ApplicationReturnStatus code) {
    switch (code) {
    case Ipopt::Solve_Succeeded: return "Solve_Succeeded";
    case Ipopt::Solved_To_Acceptable_Level: return "Solved_To_Acceptable_Level";
    case Ipopt::Infeasible_Problem_Detected: return "Infeasible_Problem_Detected";
    case Ipopt::Search_Direction_Becomes_Too_Small: return "Search_Direction_Becomes_Too_Small";
    case Ipopt::Diverging_Iterates: return "Diverging_Iterates";
    case Ipopt::User_Requested_Stop: return "User_Requested_Stop";
    case Ipopt::Feasible_Point_Found: return "Feasible_Point_Found";
    case Ipopt::Maximum_Iterations_Exceeded: return "Maximum_Iterations_Exceeded";
    case Ipopt::Restoration_Failed: return "Restoration_Failed";
    case Ipopt::Error_In_Step_Computation: return "Error_In_Step_Computation";
    case Ipopt::Maximum_CpuTime_Exceeded: return "Maximum_CpuTime_Exceeded";
    case Ipopt::Not_Enough_Degrees_Of_Freedom: return "Not_Enough_Degrees_Of_Freedom";
    case Ipopt::Invalid_Number_Detected: return "Invalid_Number_Detected";
    case Ipopt::Invalid_Problem_Definition: return "Invalid_Problem_Definition";
    case Ipopt::Invalid_Option: return "Invalid_Option";
    case Ipopt::Unrecoverable_Exception: return "Unrecoverable_Exception";
    case Ipopt::NonIpopt_Exception_Thrown: return "NonIpopt_Exception_Thrown";
    case Ipopt::Insufficient_Memory: return "Insufficient_Memory";
    case Ipopt::Internal_Error: return "Internal_Error";
    default: return std::to_string(static_cast<int>(code));
    }
}

// The one linear-solver rule, decided once per loaded library: EPCSAFT_IPOPT_LINEAR_SOLVER forces mumps or
// ma57; otherwise MA57 when Ipopt's loader opens EPCSAFT_HSL_LIBRARY (unset: libhsl.so) and finds it, else MUMPS.
inline const std::string& linear_solver() {
    static const std::string chosen = [] {
        const char* forced = std::getenv("EPCSAFT_IPOPT_LINEAR_SOLVER");
        const std::string name = forced ? forced : "";
        if (!name.empty() && name != "mumps" && name != "ma57")
            throw std::invalid_argument("EPCSAFT_IPOPT_LINEAR_SOLVER must be mumps or ma57, not " + name);
        if (name == "mumps") return name;
#if __has_include(<HSLLoader.h>)
        const char* path = std::getenv("EPCSAFT_HSL_LIBRARY");
        char message[512] = "loaded library has no MA57";
        const bool ma57 = LSL_loadHSL(path && *path ? path : nullptr, message, sizeof message) == 0
            && LSL_isMA57available() != 0;
        if (name == "ma57" && !ma57)
            throw std::invalid_argument(std::string("EPCSAFT_IPOPT_LINEAR_SOLVER=ma57 is unavailable: ") + message);
        return std::string(ma57 ? "ma57" : "mumps");
#else
        if (name == "ma57")
            throw std::invalid_argument("MA57 is unavailable: Ipopt was built without its linear-solver loader");
        return std::string("mumps");
#endif
    }();
    return chosen;
}

inline void observe_ipopt(IpoptDiagnostics& d, core::Scope* scope, const IpoptIteration& iteration) {
    ++d.callbacks;
    if (iteration.mode == int(Ipopt::RestorationPhaseMode) && (!d.last || d.last->mode != iteration.mode)) {
        ++d.restoration_entries;
        if (!d.first_restoration_iteration) d.first_restoration_iteration = iteration.iteration;
    }
    d.last = iteration;
    if (scope) scope->probe("ipopt_iteration",
        {double(iteration.mode), iteration.objective, iteration.inf_pr, iteration.inf_du, iteration.mu, iteration.step_norm, iteration.regularization,
         iteration.alpha_pr, iteration.alpha_du, double(iteration.line_search_trials)}, iteration.iteration);
}

// One checked option/lifecycle owner for exact curvature and estimation's BFGS.
inline std::pair<Ipopt::ApplicationReturnStatus, Ipopt::ApplicationReturnStatus> run_ipopt(
    const Ipopt::SmartPtr<Ipopt::TNLP>& nlp, double tolerance, int iterations,
    double seconds, bool exact, double objective_factor, std::map<std::string, double>* consumed = nullptr,
    double requested_row = 0.0, const std::function<bool(double)>& warm = {}, std::string* solver = nullptr) {
    const double dual_tolerance = tolerance * std::max(1.0, 1.0 / objective_factor);
    if (!(tolerance > 0.0) || !std::isfinite(tolerance) || !(dual_tolerance > 0.0) || !std::isfinite(dual_tolerance))
        throw std::invalid_argument("EQUILIBRIUM_ROW_RESIDUAL=" + core::optional_number(requested_row)
            + " objective_factor=" + core::optional_number(objective_factor) + " produces invalid inner target="
            + core::optional_number(tolerance) + " or dual_inf_tol=" + core::optional_number(dual_tolerance));
    const char* options_file = std::getenv("EPCSAFT_IPOPT_OPTIONS_FILE");
    // Clones share this thread's option registry and journalist (#646); an options file's output_file must not.
    thread_local const Ipopt::SmartPtr<Ipopt::IpoptApplication> base = ::IpoptApplicationFactory();
    Ipopt::SmartPtr<Ipopt::IpoptApplication> logged;
    if (auto* trace = core::current_trace(); trace && trace->ipopt_log) {
        auto& request = *trace->ipopt_log;
        using Application = Ipopt::SmartPtr<Ipopt::IpoptApplication>;
        if (!request.attempted) {
            request.attempted = true;
            auto holder = std::make_shared<Application>(::IpoptApplicationFactory());
            if ((*holder)->OpenOutputFile(request.destination, Ipopt::J_ALL)) {
                (*holder)->Jnlst()->GetJournal("OutputFile:" + request.destination)->SetPrintLevel(Ipopt::J_MAIN, Ipopt::J_ITERSUMMARY);
                request.application = std::move(holder);
            } else {
                request.error = "Ipopt log could not be opened";
                trace->incomplete = true;
            }
        }
        if (request.application) logged = *std::static_pointer_cast<Application>(request.application);
    }
    const Ipopt::SmartPtr<Ipopt::IpoptApplication> app = options_file && *options_file ? ::IpoptApplicationFactory()
        : Ipopt::IsValid(logged) ? logged->clone() : base->clone();
    const auto options = app->Options();
    const auto admit = [](bool ok) { if (!ok) throw std::invalid_argument("Ipopt option rejected"); };
    admit(options->SetIntegerValue("print_level", 0));
    admit(options->SetStringValue("sb", "yes"));  // no banner on stdout
    admit(options->SetStringValue("nlp_scaling_method", "none"));
    admit(options->SetNumericValue("obj_scaling_factor", objective_factor));
    admit(options->SetStringValue("hessian_approximation", exact ? "exact" : "limited-memory"));
    for (const auto* name : {"tol", "acceptable_tol", "constr_viol_tol"})
        admit(options->SetNumericValue(name, tolerance));
    // ALG-EQUILIBRIUM: the absolute dual stop follows the extensive objective basis.
    admit(options->SetNumericValue("dual_inf_tol", dual_tolerance));
    admit(options->SetIntegerValue("acceptable_iter", 0));
    admit(options->SetIntegerValue("max_iter", iterations));
    if (seconds > 0.0) admit(options->SetNumericValue("max_cpu_time", seconds));
    if (!exact) {
        admit(options->SetStringValue("limited_memory_update_type", "bfgs"));
        admit(options->SetNumericValue("compl_inf_tol", tolerance));
        admit(options->SetNumericValue("bound_relax_factor", 0.0));
    }
    admit(options->SetStringValue("linear_solver", linear_solver()));
#if __has_include(<HSLLoader.h>)
    admit(options->SetIntegerValue("ma57_pivot_order", 2));          // AMD, as measured in #509
    admit(options->SetStringValue("ma57_automatic_scaling", "yes")); // unscaled, a trace moiety row reads as dependent
#endif
    if (options_file && *options_file) {
        std::error_code error;
        if (!std::filesystem::is_regular_file(options_file, error) || !std::ifstream(options_file).good())
            throw std::invalid_argument(std::string("EPCSAFT_IPOPT_OPTIONS_FILE must name a readable regular file: ") + options_file);
    }
    const auto initialized = app->Initialize(options_file ? options_file : "");
    if (solver) options->GetStringValue("linear_solver", *solver, "");  // after any options-file override
    if (consumed) for (const auto* name : {"tol", "acceptable_tol", "constr_viol_tol", "dual_inf_tol"}) {
        double value = 0.0; admit(options->GetNumericValue(name, value, "")); (*consumed)[name] = value;
    }
    if (initialized == Ipopt::Solve_Succeeded && warm) {
        double cap = 0.0;
        options->GetNumericValue("warm_start_mult_init_max", cap, "");
        if (warm(options_file && *options_file ? 0.0 : cap)) {
            admit(options->SetStringValue("warm_start_init_point", "yes"));
            for (const auto* name : {"warm_start_bound_push", "warm_start_bound_frac", "warm_start_mult_bound_push",
                                    "warm_start_slack_bound_push", "warm_start_slack_bound_frac"})
                admit(options->SetNumericValue(name, policy::EQUILIBRIUM_WARM_PUSH.value));
            admit(options->SetStringValue("mu_strategy", "monotone"));
            admit(options->SetNumericValue("mu_init", policy::EQUILIBRIUM_WARM_MU.value));
        }
    }
    return {initialized, initialized == Ipopt::Solve_Succeeded ? app->OptimizeTNLP(nlp) : initialized};
}
}  // namespace epcsaft::equilibrium
