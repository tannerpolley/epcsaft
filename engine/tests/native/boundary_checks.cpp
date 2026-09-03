#include "equilibrium/problem.hpp"
#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium_support.hpp"

#include "core/numerical_policy.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

namespace {

double centered_step(double value) {
    return std::sqrt(std::numeric_limits<double>::epsilon())
        * std::max(std::abs(value), 1.0);
}

double convergence_limit(const epcsaft::core::Convergence& convergence) {
    CHECK(convergence.criterion != nullptr);
    return std::max(convergence.criterion->value, convergence.floor.value_or(0.0));
}

void check_s3() {
    epcsaft::HelmholtzModel model(model_support::neutral_binary());
    auto bubble = equilibrium_support::boundary();
    bubble.feed.values = {0.5, 0.5};
    bubble.temperature.value = 180.0;
    bubble.pressure.value = 1.0e5;
    const auto bubble_result = epcsaft::equilibrium::solve_equilibrium(model, bubble);
    CHECK(bubble_result.success);
    auto dew = bubble;
    if (!bubble_result.success) std::fprintf(stderr, "S3 bubble: %s\n", bubble_result.message.c_str());
    if (!bubble_result.success) return;
    dew.feed.values = {bubble_result.mole_fractions[2], bubble_result.mole_fractions[3]};
    dew.pressure.value = bubble_result.pressure.front();
    dew.phases[0].amount = epcsaft::equilibrium::AmountSpec::pinned(0.0);
    dew.phases[1].amount = epcsaft::equilibrium::AmountSpec::free();
    const auto dew_result = epcsaft::equilibrium::solve_equilibrium(model, dew);
    CHECK(dew_result.success);
    if (!dew_result.success) std::fprintf(stderr, "S3 dew: %s\n", dew_result.message.c_str());
    if (!dew_result.success) return;
    const double limit = std::max(convergence_limit(bubble_result.convergence),
                                  convergence_limit(dew_result.convergence));
    CHECK(std::abs(dew_result.pressure.front() / bubble_result.pressure.front() - 1.0) <= limit);
    CHECK(std::abs(dew_result.mole_fractions[0] - bubble_result.mole_fractions[0]) <= limit);
    CHECK(std::abs(dew_result.mole_fractions[1] - bubble_result.mole_fractions[1]) <= limit);
}

void check_s4() {
    epcsaft::HelmholtzModel model(model_support::neutral_ternary());
    auto bubble = equilibrium_support::flash(3);
    bubble.phases = {equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid),
                     equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor,
                                                epcsaft::equilibrium::AmountSpec::pinned(0.0))};
    bubble.temperature.value = 180.0;
    bubble.pressure = epcsaft::equilibrium::FreeValue::variable(1.0e5);
    const auto bubble_result = epcsaft::equilibrium::solve_equilibrium(model, bubble);
    CHECK(bubble_result.success);
    if (!bubble_result.success) std::fprintf(stderr, "S4 bubble: %s\n", bubble_result.message.c_str());
    if (!bubble_result.success) return;
    std::array<double, 3> gaps{};
    std::optional<epcsaft::equilibrium::EquilibriumResult> start = bubble_result;
    std::size_t index = 0;
    for (const double beta : {1.0e-2, 1.0e-4, 1.0e-6}) {
        auto specification = bubble;
        specification.phases[1].amount = epcsaft::equilibrium::AmountSpec::pinned(beta);
        specification.temperature = epcsaft::equilibrium::FreeValue::variable(180.0);
        specification.pressure = epcsaft::equilibrium::FreeValue::fixed(
            bubble_result.pressure.front());
        const auto result = epcsaft::equilibrium::solve_equilibrium(
            model, specification, start);
        CHECK(result.success && result.temperature.size() == 1);
        if (!result.success) return;
        gaps[index++] = std::abs(result.temperature.front() - bubble_result.temperature.front());
        start = result;
    }
    CHECK(gaps[1] <= gaps[0] / 10.0 && gaps[2] <= gaps[1] / 10.0);
}

void check_s5() {
    epcsaft::HelmholtzModel model(model_support::neutral_ternary());
    auto problem = equilibrium_support::flash(3);
    problem.phases = {equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid),
                      equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor,
                                                 epcsaft::equilibrium::AmountSpec::pinned(0.0))};
    problem.temperature.value = 180.0;
    problem.pressure = epcsaft::equilibrium::FreeValue::variable(1.0e5);
    const auto pressure_result = epcsaft::equilibrium::solve_equilibrium(model, problem);
    CHECK(pressure_result.success);
    if (!pressure_result.success) std::fprintf(stderr, "S5 pressure: %s\n", pressure_result.message.c_str());
    if (!pressure_result.success) return;
    problem.temperature = epcsaft::equilibrium::FreeValue::variable(
        pressure_result.temperature.front());
    problem.pressure = epcsaft::equilibrium::FreeValue::fixed(
        pressure_result.pressure.front());
    const auto temperature_result = epcsaft::equilibrium::solve_equilibrium(
        model, problem, pressure_result);
    CHECK(temperature_result.success);
    if (!temperature_result.success) return;
    CHECK(std::abs(temperature_result.temperature.front()
                   - pressure_result.temperature.front()) <= 1.0e-6);
}

void check_s9() {
    epcsaft::HelmholtzModel model(model_support::neutral_ternary());
    auto problem = equilibrium_support::flash(3);
    problem.phases = {equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor),
                      equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid)};
    problem.temperature.value = 180.0;
    const auto result = epcsaft::equilibrium::solve_equilibrium(model, problem);
    CHECK(result.success);
    if (!result.success) return;
    const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
    const auto point = result.coordinates;
    CHECK(point.size() == compiled.rows.size());
    for (std::size_t lever_index = 0; lever_index < 3; ++lever_index) {
        const double delta = lever_index == 0 ? 1.0e-3 : lever_index == 1 ? 1.0 : 1.0e-5;
        auto higher = problem;
        auto lower = problem;
        std::vector<double> lever(5, 0.0);
        lever[lever_index < 2 ? lever_index : 2] = 1.0;
        if (lever_index == 0) { higher.temperature.value += delta; lower.temperature.value -= delta; }
        if (lever_index == 1) { higher.pressure.value += delta; lower.pressure.value -= delta; }
        if (lever_index == 2) { higher.feed.values[0] += delta; lower.feed.values[0] -= delta; }
        const auto derivative = epcsaft::equilibrium::sensitivity(compiled, point, lever);
        const auto up = epcsaft::equilibrium::solve_equilibrium(model, higher, result);
        const auto down = epcsaft::equilibrium::solve_equilibrium(model, lower, result);
        CHECK(up.success && down.success);
        if (!up.success || !down.success) continue;
        for (std::size_t column = 0; column < point.size(); ++column) {
            const double centered = (up.coordinates[column] - down.coordinates[column]) / (2.0 * delta);
            CHECK(std::abs(derivative[column] - centered) <= 1.0e-6
                  * std::max(1.0, std::abs(centered)));
        }
    }
}

void check_s10() {
    epcsaft::HelmholtzModel model(model_support::neutral_binary());
    const auto compiled = epcsaft::equilibrium::compile_problem(
        model, equilibrium_support::boundary());
    const auto point = compiled.generic_point();
    const auto jacobian = compiled.jacobian(point);
    const auto values = compiled.evaluate(point).values;
    CHECK(jacobian.size() == values.size() * point.size());
    for (std::size_t column = 0; column < point.size(); ++column) {
        auto plus = point; auto minus = point;
        const double delta = centered_step(point[column]);
        plus[column] += delta; minus[column] -= delta;
        const auto up = compiled.evaluate(plus).values;
        const auto down = compiled.evaluate(minus).values;
        for (std::size_t row = 0; row < values.size(); ++row)
            CHECK(std::abs(jacobian[row * point.size() + column]
                          - (up[row] - down[row]) / (2.0 * delta))
                  <= 1.0e-6 * std::max(1.0, std::abs(up[row] - down[row]) / (2.0 * delta)));
    }
}

void check_s14() {
    epcsaft::HelmholtzModel model(model_support::ethane());
    auto problem = equilibrium_support::flash(1);
    problem.phases = {equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid),
                      equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor,
                                                 epcsaft::equilibrium::AmountSpec::pinned(0.0))};
    problem.temperature.value = 250.0;
    problem.pressure = epcsaft::equilibrium::FreeValue::variable(1.3e6);
    const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
    CHECK(compiled.rank == 1 && compiled.generic_nullity == 0);
    CHECK(compiled.charge_row_count == 0);
    CHECK(compiled.rows_minus_unknowns == 0);
    for (const auto& row : compiled.rows) CHECK(row.kind != epcsaft::equilibrium::RowKind::Electro);
    for (const auto& row : compiled.rows) std::printf("S14 row %s\n", row.eqid);
    const auto vapor_pinned = epcsaft::equilibrium::solve_equilibrium(model, problem);
    CHECK(vapor_pinned.success);
    if (!vapor_pinned.success) return;
    problem.phases[0].amount = epcsaft::equilibrium::AmountSpec::pinned(0.0);
    problem.phases[1].amount = epcsaft::equilibrium::AmountSpec::free();
    const auto liquid_pinned = epcsaft::equilibrium::solve_equilibrium(model, problem, vapor_pinned);
    CHECK(liquid_pinned.success);
    if (!liquid_pinned.success) return;
    const double limit = std::max(convergence_limit(vapor_pinned.convergence),
                                  convergence_limit(liquid_pinned.convergence));
    CHECK(std::abs(vapor_pinned.pressure.front() / liquid_pinned.pressure.front() - 1.0) <= limit);
}

}  // namespace

int main() {
    equilibrium_support::run_check("S3", [] { check_s3(); });
    equilibrium_support::run_check("S4", [] { check_s4(); });
    equilibrium_support::run_check("S5", [] { check_s5(); });
    equilibrium_support::run_check("S9", [] { check_s9(); });
    equilibrium_support::run_check("S10", [] { check_s10(); });
    equilibrium_support::run_check("S14", [] { check_s14(); });
    return equilibrium_support::failures == 0 ? 0 : 1;
}
