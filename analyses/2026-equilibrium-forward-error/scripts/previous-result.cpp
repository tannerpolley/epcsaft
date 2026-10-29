#include "equilibrium_support.hpp"
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace epcsaft;
using namespace epcsaft::equilibrium;

// The delivered-review failing call: identical declaration, then reuse its refined result.
int main() {
    HelmholtzModel model(properties_support::neutral_binary());
    auto problem = equilibrium_support::bubble_at(265.95);
    problem.pressure = FreeValue::variable(6.75e6);
    EquilibriumResult original;
    std::cout << std::setprecision(20)
        << "capture,P,rhoL,rhoV,radiusL,radiusV,pressure_radius,bits,iterations,evaluations,hessians,corrections,activations,decisions,first_failure\n";
    for (bool detailed : {false, true, false}) {
        std::optional<core::Capture> capture;
        if (detailed) capture.emplace(true);
        const auto first = solve_equilibrium(model, problem);
        const auto result = solve_equilibrium(model, problem, first);
        if (!result.success || !result.numerical_estimates || !result.linear_solve_diagnostics)
            throw std::runtime_error("previous-result endpoint unavailable");
        const auto& e = *result.numerical_estimates;
        constexpr long double rhoL = 8676.4528556443341908204224L, rhoV = 8665.1463344039600528343775L;
        if (std::abs(result.molar_densities[1] / rhoV - 1) > 1e-7L
            || std::abs(std::log(result.molar_densities[0] / rhoL)) > e.densities[0].chart
            || std::abs(std::log(result.molar_densities[1] / rhoV)) > e.densities[1].chart)
            throw std::runtime_error("previous-result density accuracy/estimate falsified");
        if (!original.coordinates.empty() && (original.coordinates != result.coordinates
            || original.evaluations != result.evaluations || original.iterations != result.iterations
            || original.hessian_calls != result.hessian_calls
            || original.numerical_estimates->physical_radius != e.physical_radius))
            throw std::runtime_error("capture changed state, estimate or work");
        original = result;
        std::size_t activations = 0, decisions = 0;
        if (capture) for (const auto& event : capture->trace().events) {
            activations += std::string_view(event.name) == "mixture_precision_activate";
            decisions += std::string_view(event.name) == "mixture_precision_decision";
        }
        if (capture && (activations != 2 || decisions != 2))
            throw std::runtime_error("both ordinary calls must select and activate precision");
        std::cout << detailed << ',' << *result.pressure;
        for (double rho : result.molar_densities) std::cout << ',' << rho;
        std::cout << ',' << e.densities[0].chart << ',' << e.densities[1].chart << ',' << e.pressure.absolute
            << ',' << result.linear_solve_diagnostics->precision_bits << ',' << result.iterations
            << ',' << result.evaluations << ',' << result.hessian_calls
            << ',' << result.linear_solve_diagnostics->refinement_iterations << ',' << activations << ',' << decisions
            << ',' << std::quoted(result.first_numerical_failure ? result.first_numerical_failure->cause : "") << '\n';
    }
}
