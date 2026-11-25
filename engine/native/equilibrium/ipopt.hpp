#include "core/issue235_meter.hpp"
#pragma once
#include "core/convergence.hpp"
#include "equilibrium/problem.hpp"
#include "core/trace.hpp"
#include <IpLibraryLoader.hpp>
#include <IpLinearSolvers.h>
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
        if (IpoptGetAvailableLinearSolvers(1) & IPOPTLINEARSOLVER_MA57) return std::string("ma57");
        const char* path = std::getenv("EPCSAFT_HSL_LIBRARY");
        try {
            Ipopt::LibraryLoader loader(path && *path ? path : "libhsl.so");
            loader.loadSymbol("ma57ad");
            return std::string("ma57");
        } catch (const Ipopt::DYNAMIC_LIBRARY_FAILURE& error) {
            if (name == "ma57") throw std::invalid_argument(error.Message());
            return std::string("mumps");
        }
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
    if (auto* log = ipopt_iteration_log())
        log->push_back({(log->empty() ? 0 : log->back().run) + (iteration.iteration == 0), iteration});
    if (scope) scope->probe("ipopt_iteration",
        {double(iteration.mode), iteration.objective, iteration.inf_pr, iteration.inf_du, iteration.mu, iteration.step_norm, iteration.regularization,
         iteration.alpha_pr, iteration.alpha_du, double(iteration.line_search_trials)}, iteration.iteration);
}

// One checked option/lifecycle owner for exact curvature and estimation's BFGS.
inline std::pair<Ipopt::ApplicationReturnStatus, Ipopt::ApplicationReturnStatus> run_ipopt(
    const Ipopt::SmartPtr<Ipopt::TNLP>& nlp, double tolerance, int iterations,
    double seconds, bool exact, double objective_factor, std::map<std::string, double>* consumed = nullptr,
    double requested_row = 0.0, const std::function<bool(double)>& warm = {}, std::string* solver = nullptr, bool amount_limited = false, bool feasibility = false, bool free_pressure = false, bool nonreacting_multiphase = false) {
    const double dual_tolerance = tolerance * std::max(1.0, 1.0 / objective_factor);
    if (!(tolerance > 0.0) || !std::isfinite(tolerance) || !(dual_tolerance > 0.0) || !std::isfinite(dual_tolerance))
        throw std::invalid_argument("EQUILIBRIUM_ROW_RESIDUAL=" + core::optional_number(requested_row)
            + " objective_factor=" + core::optional_number(objective_factor) + " produces invalid inner target="
            + core::optional_number(tolerance) + " or dual_inf_tol=" + core::optional_number(dual_tolerance));
    ep235::Scope lifecycle_meter("ipopt_setup");
    core::Scope lifecycle("ipopt_application");
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
            if ((*holder)->OpenOutputFile(request.destination, Ipopt::EJournalLevel(request.print_level))) {
                (*holder)->Jnlst()->GetJournal("OutputFile:" + request.destination)->SetPrintLevel(Ipopt::J_MAIN,
                    Ipopt::EJournalLevel(std::min(request.print_level, int(Ipopt::J_ITERSUMMARY))));
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
    const auto admit = [](bool ok, const char* name) { if (!ok) throw std::invalid_argument(std::string("Ipopt option rejected: ") + name); };
    admit(options->SetIntegerValue("print_level", 0), "print_level");
    admit(options->SetStringValue("sb", "yes"), "sb");  // no banner on stdout
    admit(options->SetStringValue("nlp_scaling_method", "none"), "nlp_scaling_method");
    admit(options->SetNumericValue("obj_scaling_factor", objective_factor), "obj_scaling_factor");
    admit(options->SetStringValue("hessian_approximation", exact ? "exact" : "limited-memory"), "hessian_approximation");
    for (const auto* name : {"tol", "acceptable_tol", "constr_viol_tol"})
        admit(options->SetNumericValue(name, tolerance), name);
    // ALG-EQUILIBRIUM: the absolute dual stop follows the extensive objective basis.
    admit(options->SetNumericValue("dual_inf_tol", dual_tolerance), "dual_inf_tol");
    admit(options->SetIntegerValue("acceptable_iter", 0), "acceptable_iter");
    admit(options->SetIntegerValue("max_iter", iterations), "max_iter");
    if (seconds > 0.0) admit(options->SetNumericValue("max_cpu_time", seconds), "max_cpu_time");
    // ALG-EQUILIBRIUM: unrelaxed variable bounds; a relaxed zero amount lets a vanished phase lower the objective.
    admit(options->SetNumericValue("bound_relax_factor", 0.0), "bound_relax_factor");
    if (amount_limited && (feasibility ? !(free_pressure && linear_solver() == "ma57") : !nonreacting_multiphase)) {
        admit(options->SetNumericValue("bound_relax_factor", feasibility ? 1e-8 : 1e-12), "bound_relax_factor");
    }
    if (!exact) {
        admit(options->SetStringValue("limited_memory_update_type", "bfgs"), "limited_memory_update_type");
        admit(options->SetNumericValue("compl_inf_tol", tolerance), "compl_inf_tol");
    }
    admit(options->SetStringValue("linear_solver", linear_solver()), "linear_solver");
    if (linear_solver() == "ma57") {
        if (const char* path = std::getenv("EPCSAFT_HSL_LIBRARY"); path && *path)
            admit(options->SetStringValue("hsllib", path), "hsllib");
        admit(options->SetIntegerValue("ma57_pivot_order", 2), "ma57_pivot_order");          // AMD, as measured in #509
        admit(options->SetStringValue("ma57_automatic_scaling", "yes"), "ma57_automatic_scaling"); // unscaled, a trace moiety row reads as dependent
    }
    if (options_file && *options_file) {
        std::error_code error;
        if (!std::filesystem::is_regular_file(options_file, error) || !std::ifstream(options_file).good())
            throw std::invalid_argument(std::string("EPCSAFT_IPOPT_OPTIONS_FILE must name a readable regular file: ") + options_file);
    }
    if (feasibility) {
        if (linear_solver() == "ma57") admit(options->SetStringValue("mu_strategy", free_pressure ? "adaptive" : "monotone"), "mu_strategy");
        else admit(options->SetStringValue("recalc_y", "yes"), "recalc_y");
    }
    // Lock both triggers before Initialize reads a diagnostic file: rollback invalidates the saved amount base.
    if (amount_limited) for (const auto* name : {"watchdog_shortened_iter_trigger", "resto.watchdog_shortened_iter_trigger"})
        admit(options->SetIntegerValue(name, 0, false), name);
    const auto initialized = app->Initialize(options_file ? options_file : "");
    if (solver) options->GetStringValue("linear_solver", *solver, "");  // after any options-file override
    if (consumed) for (const auto* name : {"tol", "acceptable_tol", "constr_viol_tol", "dual_inf_tol"}) {
        double value = 0.0; admit(options->GetNumericValue(name, value, ""), name); (*consumed)[name] = value;
    }
    if (consumed) options->GetNumericValue("bound_push", (*consumed)["bound_push"], "");  // populated when unset; false then
    if (initialized == Ipopt::Solve_Succeeded && warm) {
        double cap = 0.0;
        options->GetNumericValue("warm_start_mult_init_max", cap, "");
        if (warm(options_file && *options_file ? 0.0 : cap)) {
            admit(options->SetStringValue("warm_start_init_point", "yes"), "warm_start_init_point");
            for (const auto* name : {"warm_start_bound_push", "warm_start_bound_frac", "warm_start_mult_bound_push",
                                    "warm_start_slack_bound_push", "warm_start_slack_bound_frac"})
                admit(options->SetNumericValue(name, policy::EQUILIBRIUM_WARM_PUSH.value), name);
            admit(options->SetStringValue("mu_strategy", "monotone"), "mu_strategy");
            admit(options->SetNumericValue("mu_init", policy::EQUILIBRIUM_WARM_MU.value), "mu_init");
        }
    }
    if (initialized != Ipopt::Solve_Succeeded) { lifecycle_meter.key.outcome = static_cast<int>(initialized); return {initialized, initialized}; }
    ep235::Scope optimizer_meter("ipopt_optimize");
    core::Scope optimize("ipopt_optimize");
    const auto status = app->OptimizeTNLP(nlp);
    optimizer_meter.key.outcome = static_cast<int>(status);
    return {initialized, status};
}
}  // namespace epcsaft::equilibrium
