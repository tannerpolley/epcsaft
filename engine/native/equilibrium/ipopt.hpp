#pragma once
#include <IpIpoptApplication.hpp>
#include <IpTNLP.hpp>
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
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

// One checked option/lifecycle owner for exact curvature and estimation's BFGS.
inline std::pair<Ipopt::ApplicationReturnStatus, Ipopt::ApplicationReturnStatus> run_ipopt(
    const Ipopt::SmartPtr<Ipopt::TNLP>& nlp, double tolerance, int iterations,
    double seconds, bool exact, double objective_factor) {
    auto app = Ipopt::SmartPtr<Ipopt::IpoptApplication>(::IpoptApplicationFactory());
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
    const double dual_tolerance = tolerance * std::max(1.0, 1.0 / objective_factor);
    admit(options->SetNumericValue("dual_inf_tol", dual_tolerance));
    admit(options->SetIntegerValue("acceptable_iter", 0));
    admit(options->SetIntegerValue("max_iter", iterations));
    if (seconds > 0.0) admit(options->SetNumericValue("max_cpu_time", seconds));
    if (!exact) {
        admit(options->SetStringValue("limited_memory_update_type", "bfgs"));
        admit(options->SetNumericValue("compl_inf_tol", tolerance));
        admit(options->SetNumericValue("bound_relax_factor", 0.0));
    }
    const char* options_file = std::getenv("EPCSAFT_IPOPT_OPTIONS_FILE");
    if (options_file && *options_file) {
        std::error_code error;
        if (!std::filesystem::is_regular_file(options_file, error) || !std::ifstream(options_file).good())
            throw std::invalid_argument(std::string("EPCSAFT_IPOPT_OPTIONS_FILE must name a readable regular file: ") + options_file);
    }
    const auto initialized = app->Initialize(options_file ? options_file : "");
    return {initialized, initialized == Ipopt::Solve_Succeeded ? app->OptimizeTNLP(nlp) : initialized};
}
}  // namespace epcsaft::equilibrium
