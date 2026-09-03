#include "equilibrium/problem.hpp"
#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium_support.hpp"

#include "core/numerical_policy.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

#define CHECK(condition) do { if (!(condition)) { std::fprintf(stderr, "%s:%d: failed: %s\n", __FILE__, __LINE__, #condition); std::abort(); } } while (false)

namespace {

double centered_step(double value) {
    return std::sqrt(std::numeric_limits<double>::epsilon())
        * std::max(std::abs(value), 1.0);
}

void check_s3() {
    epcsaft::HelmholtzModel model(model_support::associating_binary());
    auto bubble = equilibrium_support::boundary();
    bubble.feed.values = {0.5, 0.5};
    bubble.temperature.value = 280.0;
    bubble.pressure.value = 2500.0;
    const auto bubble_result = epcsaft::equilibrium::solve_equilibrium(model, bubble);
    CHECK(bubble_result.success);
    auto dew = bubble;
    if (bubble_result.success) {
        dew.feed.values = {bubble_result.mole_fractions[2], bubble_result.mole_fractions[3]};
    }
    dew.phases[0].amount = epcsaft::equilibrium::AmountSpec::pinned(0.0);
    dew.phases[1].amount = epcsaft::equilibrium::AmountSpec::free();
    if (bubble_result.success) {
        CHECK(bubble_result.pressure.size() == 1);
        const auto dew_result = epcsaft::equilibrium::solve_equilibrium(model, dew, bubble_result);
        CHECK(dew_result.success);
        {
            CHECK(dew_result.has_state());
            CHECK(dew_result.convergence.criterion != nullptr);
            CHECK(std::isfinite(dew_result.convergence.residual));
            const double limit = 1.0e-6 * std::max(1.0, std::abs(bubble_result.pressure.front()));
            CHECK(std::abs(dew_result.pressure.front() - bubble_result.pressure.front()) < limit);
            CHECK(std::abs(dew_result.mole_fractions[2] - bubble_result.mole_fractions[0]) < limit);
            CHECK(std::abs(dew_result.mole_fractions[3] - bubble_result.mole_fractions[1]) < limit);
        }
    }
}

void check_s4() {
    epcsaft::HelmholtzModel model(model_support::neutral_binary());
    auto bubble = equilibrium_support::boundary();
    bubble.temperature.value = 280.0;
    bubble.pressure.value = 2500.0;
    const auto bubble_result = epcsaft::equilibrium::solve_equilibrium(model, bubble);
    CHECK(bubble_result.success);
    std::array<double, 3> gaps{};
    std::optional<epcsaft::equilibrium::EquilibriumResult> start = bubble_result;
    std::size_t index = 0;
    for (const double beta : {1.0e-2, 1.0e-4, 1.0e-6}) {
        auto specification = bubble;
        specification.phases[1].amount = epcsaft::equilibrium::AmountSpec::pinned(beta);
        specification.temperature = epcsaft::equilibrium::FreeValue::variable(280.0);
        specification.pressure = epcsaft::equilibrium::FreeValue::fixed(
            bubble_result.pressure.front());
        const auto result = epcsaft::equilibrium::solve_equilibrium(
            model, specification, start);
        CHECK(result.success && result.temperature.size() == 1);
        gaps[index++] = std::abs(result.temperature.front() - bubble_result.temperature.front());
        start = result;
    }
    CHECK(gaps[1] <= gaps[0] / 10.0 && gaps[2] <= gaps[1] / 10.0);
}

void check_s5() {
    epcsaft::HelmholtzModel model(model_support::neutral_binary());
    auto problem = equilibrium_support::boundary();
    problem.temperature.value = 280.0;
    problem.pressure.value = 2500.0;
    const auto pressure_result = epcsaft::equilibrium::solve_equilibrium(model, problem);
    CHECK(pressure_result.success);
    problem.temperature = epcsaft::equilibrium::FreeValue::variable(
        pressure_result.temperature.front());
    problem.pressure = epcsaft::equilibrium::FreeValue::fixed(
        pressure_result.pressure.front());
    const auto temperature_result = epcsaft::equilibrium::solve_equilibrium(
        model, problem, pressure_result);
    CHECK(temperature_result.success);
    CHECK(std::abs(temperature_result.temperature.front()
                   - pressure_result.temperature.front()) < 1.0e-6);
}

void check_s9() {
    epcsaft::HelmholtzModel model(model_support::neutral_binary());
    auto problem = equilibrium_support::flash(2);
    problem.temperature.value = 280.0;
    const auto result = epcsaft::equilibrium::solve_equilibrium(model, problem);
    CHECK(result.success);
    const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
    const auto point = result.coordinates;
    CHECK(point.size() == compiled.rows.size());
    const double delta = 1.0e-3;
    auto higher = problem;
    higher.temperature.value += delta;
    auto lower = problem;
    lower.temperature.value -= delta;
    const auto derivative = epcsaft::equilibrium::sensitivity(compiled, point, {1.0, 0.0});
    const auto higher_result = epcsaft::equilibrium::solve_equilibrium(model, higher);
    const auto lower_result = epcsaft::equilibrium::solve_equilibrium(model, lower);
    CHECK(higher_result.success && lower_result.success);
    const auto higher_point = higher_result.coordinates;
    const auto lower_point = lower_result.coordinates;
    CHECK(derivative.size() == point.size());
    for (std::size_t column = 0; column < point.size(); ++column) {
        const double centered = (higher_point[column] - lower_point[column]) / (2.0 * delta);
        CHECK(std::abs(derivative[column] - centered) < 1.0e-6
              * std::max(1.0, std::abs(centered)));
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
    const auto problem = equilibrium_support::flash(1);
    const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
    CHECK(compiled.rank == 1);
    CHECK(compiled.charge_row_count == 0);
    CHECK(compiled.rows_minus_unknowns == 0);
    for (const auto& row : compiled.rows) CHECK(row.kind != epcsaft::equilibrium::RowKind::Electro);
}

}  // namespace

int main() {
    check_s3();
    check_s4();
    check_s5();
    check_s9();
    check_s10();
    check_s14();
    std::printf("boundary checks passed\n");
}
