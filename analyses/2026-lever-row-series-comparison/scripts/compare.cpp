#include "equilibrium_support.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
using namespace epcsaft;
using namespace epcsaft::equilibrium;
volatile double sink;
template<class Function> double cost(Function f) {
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < 30; ++i) { auto values = f(); sink = values.back(); }
    return std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - start).count() / 30;
}
void compare(const char* name, const HelmholtzModel& model, const EquilibriumProblem& problem) {
    const auto compiled = compile_problem(model, problem);
    const auto central = std::string(name) == "ionic-packing-diameter" ? solve_equilibrium(model, problem) : EquilibriumResult{};
    if (std::string(name) == "ionic-packing-diameter" && !central.success)
        throw std::runtime_error("ionic equilibrium not accepted");
    const auto point = central.success ? central.coordinates : compiled.generic_point();
    for (std::size_t i = 0; i < point.size(); ++i) std::cout << "center," << name << "," << i << ",0," << point[i] << ",0,0,0\n";
    const auto columns = compiled.lever_jacobian(point);
    const std::size_t width = columns.size() / compiled.rows.size();
    for (std::size_t i = 0; i < columns.size(); ++i)
        std::cout << "value," << name << ',' << i / width << ',' << i % width << ',' << columns[i] << '\n';
    (void)cost([&]{return compiled.lever_jacobian(point);});
    for (int repetition = 0; repetition < 7; ++repetition)
        std::cout << "cost," << name << ',' << repetition << ','
                  << cost([&]{return compiled.lever_jacobian(point);}) << '\n';
}
int main() {
    std::cout << std::setprecision(17);
    auto neutral = equilibrium_support::flash(3);
    HelmholtzModel neutral_model(properties_support::neutral_ternary(), {{ParameterFamily::Kij,{1,3}},{ParameterFamily::SegmentCount,{0}}});
    compare("neutral-fixed", neutral_model, neutral);
    neutral.phases.push_back(equilibrium_support::phase("V", PhaseKind::Vapor, AmountSpec::pinned(0.2)));
    neutral.phases.push_back(equilibrium_support::phase("L2", PhaseKind::Liquid, AmountSpec::pinned(0.3)));
    neutral.temperature = FreeValue::variable(300);
    neutral.pressure = FreeValue::variable(1e5);
    neutral.feed.values = {0.6, 0.4, 1.0};
    compare("neutral-free-TP-pinned-scaled", neutral_model, neutral);
    auto ionic = equilibrium_support::flash(3);
    ionic.feed.values = {0.98, 0.01, 0.01};
    HelmholtzModel ionic_model(properties_support::aqueous_electrolyte(), {{ParameterFamily::PackingDiameter,{1}}});
    compare("ionic-packing-diameter", ionic_model, ionic);
}
