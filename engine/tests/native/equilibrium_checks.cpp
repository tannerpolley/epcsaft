#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium/stoichiometry.hpp"
#include "equilibrium_support.hpp"
#include "core/numerical_policy.hpp"
#include "core/physical_constants.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

#define CHECK(condition) do { if (!(condition)) { std::fprintf(stderr, "%s:%d: failed: %s\n", __FILE__, __LINE__, #condition); std::abort(); } } while (false)

namespace {

double convergence_limit(const epcsaft::core::Convergence& convergence) {
    CHECK(convergence.criterion != nullptr);
    return std::max(convergence.criterion->value, convergence.floor.value_or(0.0));
}

bool close(double left, double right) {
    return std::abs(left - right) <= 1.0e-6 * std::max(1.0, std::abs(right));
}

double centered_step(double value) {
    return std::sqrt(std::numeric_limits<double>::epsilon())
        * std::max(std::abs(value), 1.0);
}

void check_finite_difference(const epcsaft::equilibrium::CompiledSystem& compiled,
                             const std::vector<double>& point) {
    const auto analytic = compiled.jacobian(point);
    const auto base = compiled.evaluate(point).values;
    CHECK(analytic.size() == base.size() * point.size());
    for (std::size_t column = 0; column < point.size(); ++column) {
        auto plus = point;
        auto minus = point;
        const double delta = centered_step(point[column]);
        plus[column] += delta;
        minus[column] -= delta;
        const auto up = compiled.evaluate(plus).values;
        const auto down = compiled.evaluate(minus).values;
        for (std::size_t row = 0; row < base.size(); ++row)
            CHECK(close(analytic[row * point.size() + column],
                        (up[row] - down[row]) / (2.0 * delta)));
    }
}

epcsaft::equilibrium::EquilibriumResult single_phase() {
    epcsaft::HelmholtzModel model(model_support::neutral_binary());
    return epcsaft::equilibrium::solve_equilibrium(model, equilibrium_support::flash(2));
}

void check_s1() {
    const auto result = single_phase();
    CHECK(result.success);
    const double limit = convergence_limit(result.convergence);
    CHECK(result.residuals.size() == result.rows.size());
    for (const double residual : result.residuals) CHECK(std::abs(residual) <= limit);
}

void check_s2() {
    epcsaft::HelmholtzModel model(model_support::associating_binary());
    auto binary = equilibrium_support::flash(2);
    binary.phases = {equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor),
                     equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid)};
    binary.phases[0].composition_guess = {0.95, 0.05};
    binary.phases[1].composition_guess = {0.05, 0.95};
    binary.temperature.value = 280.0;
    binary.pressure.value = 2500.0;
    const auto result = epcsaft::equilibrium::solve_equilibrium(model, binary);
    CHECK(result.success && result.reduced_hessian_minimum > 0.0);
    const auto reactive = epcsaft::equilibrium::solve_equilibrium(
        epcsaft::HelmholtzModel(equilibrium_support::ascani_reactive_ternary()),
        equilibrium_support::ascani_reactive_problem());
    CHECK(reactive.success && reactive.reduced_hessian_minimum > 0.0);
}

void check_s6() {
    epcsaft::HelmholtzModel model(model_support::associating_binary());
    auto one = equilibrium_support::flash(2);
    auto two = one;
    two.phases = {
        equilibrium_support::phase("L1", epcsaft::equilibrium::PhaseKind::Liquid),
        equilibrium_support::phase("L2", epcsaft::equilibrium::PhaseKind::Liquid),
    };
    const auto one_system = epcsaft::equilibrium::compile_problem(model, one);
    const auto two_system = epcsaft::equilibrium::compile_problem(model, two);
    const auto one_value = one_system.evaluate(one_system.generic_point()).objective;
    const auto two_value = two_system.evaluate(two_system.generic_point()).objective;
    CHECK(std::isfinite(one_value) && std::isfinite(two_value));
    CHECK(two_value <= one_value + 1.0e6);
}

void check_s7() {
    epcsaft::HelmholtzModel model(model_support::neutral_binary());
    auto left = equilibrium_support::flash(2);
    auto right = left;
    left.feed.values = {0.2, 0.8};
    right.feed.values = {0.4, 1.6};
    const auto a = epcsaft::equilibrium::compile_problem(model, left);
    const auto b = epcsaft::equilibrium::compile_problem(model, right);
    for (std::size_t i = 0; i < a.feed_composition.size(); ++i)
        CHECK(std::abs(a.feed_composition[i] - b.feed_composition[i]) < 1.0e-12);
}

void check_s8() {
    epcsaft::HelmholtzModel model(model_support::neutral_binary());
    auto small = equilibrium_support::flash(2);
    auto large = small;
    large.feed.values = {5.0, 5.0};
    const auto first = epcsaft::equilibrium::solve_equilibrium(model, small);
    const auto second = epcsaft::equilibrium::solve_equilibrium(model, large);
    CHECK(first.success && second.success);
    CHECK(first.temperature == second.temperature);
    CHECK(first.pressure == second.pressure);
    CHECK(first.amounts.size() == second.amounts.size());
    for (std::size_t i = 0; i < first.amounts.size(); ++i)
        CHECK(std::abs(second.amounts[i] / first.amounts[i] - 10.0) < 1.0e-7);
}

void check_s10() {
    epcsaft::HelmholtzModel model(model_support::neutral_binary());
    const auto flash = epcsaft::equilibrium::compile_problem(
        model, equilibrium_support::flash(2));
    const auto flash_point = flash.generic_point();
    check_finite_difference(flash, flash_point);
    std::vector<double> primal(flash_point.begin(),
                               flash_point.begin() + static_cast<std::ptrdiff_t>(flash.primal_size));
    epcsaft::equilibrium::GibbsMinimizer nlp(flash, primal);
    Ipopt::Index n = 0, constraints = 0, jacobian_entries = 0, hessian_entries = 0;
    Ipopt::TNLP::IndexStyleEnum style = Ipopt::TNLP::C_STYLE;
    CHECK(nlp.get_nlp_info(n, constraints, jacobian_entries, hessian_entries, style));
    std::vector<Ipopt::Index> rows(static_cast<std::size_t>(hessian_entries));
    std::vector<Ipopt::Index> columns(static_cast<std::size_t>(hessian_entries));
    std::vector<Ipopt::Number> lambda(static_cast<std::size_t>(constraints), 0.0);
    CHECK(nlp.eval_h(n, primal.data(), true, 1.0, constraints, lambda.data(),
                     true, hessian_entries, rows.data(), columns.data(), nullptr));
    std::vector<Ipopt::Number> hessian(static_cast<std::size_t>(hessian_entries), 0.0);
    CHECK(nlp.eval_h(n, primal.data(), false, 1.0, constraints, lambda.data(),
                     false, hessian_entries, rows.data(), columns.data(), hessian.data()));
    for (std::size_t column = 0; column < primal.size(); ++column) {
        auto plus = primal;
        auto minus = primal;
        const double step = centered_step(primal[column]);
        plus[column] += step;
        minus[column] -= step;
        std::vector<Ipopt::Number> grad_plus(static_cast<std::size_t>(n));
        std::vector<Ipopt::Number> grad_minus(static_cast<std::size_t>(n));
        CHECK(nlp.eval_grad_f(n, plus.data(), true, grad_plus.data()));
        CHECK(nlp.eval_grad_f(n, minus.data(), true, grad_minus.data()));
        for (std::size_t row = column; row < primal.size(); ++row) {
            const double centered = (grad_plus[row] - grad_minus[row]) / (2.0 * step);
            const std::size_t index = row * (row + 1) / 2 + column;
            CHECK(close(hessian[index], centered));
        }
    }
}

void check_s11() {
    const std::vector<epcsaft::equilibrium::ReactionSpec> reactions = {
        {"isomerisation", {-1.0, 1.0, 0.0}, 2.0, std::nullopt}
    };
    const auto standard = epcsaft::equilibrium::standard_state(300.0, 3, reactions, {});
    const auto affinity = epcsaft::equilibrium::reaction_affinities(reactions, standard);
    CHECK(affinity.size() == 1);
    CHECK(std::isfinite(affinity.front()));
    CHECK(std::abs(affinity.front() + epcsaft::GAS_CONSTANT_J_PER_MOL_K * 300.0
                   * std::log(2.0)) < 1.0e-10);
    epcsaft::HelmholtzModel model(equilibrium_support::ascani_reactive_ternary());
    auto problem = equilibrium_support::flash(3);
    problem.phases.front().kind = epcsaft::equilibrium::PhaseKind::Vapor;
    problem.reactions = reactions;
    double previous = std::numeric_limits<double>::infinity();
    for (const double pressure : {1.0e4, 1.0e3, 1.0e2}) {
        problem.pressure.value = pressure;
        const auto result = epcsaft::equilibrium::solve_equilibrium(model, problem);
        CHECK(result.success && result.mole_fractions.size() == 3);
        const double gap = std::abs(std::log(result.mole_fractions[1] / result.mole_fractions[0])
            + affinity.front() / (epcsaft::GAS_CONSTANT_J_PER_MOL_K * result.temperature.front()));
        CHECK(gap <= previous);
        previous = gap;
        if (pressure == 1.0e2) CHECK(gap <= 1.0e-3);
    }
}

void check_s13() {
    epcsaft::HelmholtzModel model(model_support::neutral_binary());
    const auto problem = equilibrium_support::flash(2);
    const auto result = epcsaft::equilibrium::solve_equilibrium(model, problem);
    CHECK(result.success && result.phase_count == 1);
    const double limit = convergence_limit(result.convergence);
    CHECK(result.mole_fractions.size() == 2);
    CHECK(std::abs(result.mole_fractions[0] - problem.feed.values[0]) <= limit);
    CHECK(std::abs(result.mole_fractions[1] - problem.feed.values[1]) <= limit);
    CHECK(result.molar_densities.size() == 1 && result.molar_densities.front() > 0.0);
    const auto root = epcsaft::solve_density(
        model, 300.0, 1.0e5, problem.feed.values, epcsaft::DensityBranch::Liquid, std::nullopt);
    const double density_limit = std::max(limit, convergence_limit(root.convergence));
    CHECK(std::abs(1.0 / result.molar_densities.front() - 1.0 / root.molar_density)
          <= density_limit);
}

void check_s15() {
    epcsaft::HelmholtzModel model(model_support::neutral_binary());
    auto problem = equilibrium_support::flash(2);
    problem.phases = {
        equilibrium_support::phase("L1", epcsaft::equilibrium::PhaseKind::Liquid),
        equilibrium_support::phase("L2", epcsaft::equilibrium::PhaseKind::Liquid),
    };
    const auto result = epcsaft::equilibrium::solve_equilibrium(model, problem);
    if (result.success) {
        CHECK(result.mole_fractions.size() == 4);
        CHECK(std::abs(result.mole_fractions[0] - result.mole_fractions[2]) >
              epcsaft::policy::STABILITY_COMPOSITION_FLOOR.value);
    } else {
        CHECK(!result.has_state());
        CHECK(!result.message.empty());
    }
}

void check_s16() {
    epcsaft::HelmholtzModel model(model_support::neutral_binary());
    auto first = equilibrium_support::flash(2);
    auto second = first;
    first.feed.values = {0.25, 0.75};
    second.feed.values = {0.50, 1.50};
    const auto left = epcsaft::equilibrium::compile_problem(model, first);
    const auto right = epcsaft::equilibrium::compile_problem(model, second);
    CHECK(left.rank == right.rank && left.balance.size() == right.balance.size());
    CHECK(epcsaft::equilibrium::vector_norm(right.n_feed)
          > epcsaft::equilibrium::vector_norm(left.n_feed));
}

}  // namespace

int main() {
    check_s1();
    check_s2();
    check_s6();
    check_s7();
    check_s8();
    check_s10();
    check_s11();
    check_s13();
    check_s15();
    check_s16();
    std::printf("equilibrium checks passed\n");
}
