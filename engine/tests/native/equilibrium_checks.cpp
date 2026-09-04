#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium/stoichiometry.hpp"
#include "equilibrium_support.hpp"
#include "core/numerical_policy.hpp"
#include "core/physical_constants.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

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

void verify_s1(const epcsaft::HelmholtzModel& model,
               const epcsaft::equilibrium::EquilibriumProblem& problem,
               const epcsaft::equilibrium::EquilibriumResult& result) {
    CHECK(result.success);
    if (!result.success) return;
    const double limit = convergence_limit(result.convergence);
    CHECK(result.residuals.size() == result.rows.size());
    for (const double residual : result.residuals) CHECK(std::abs(residual) <= limit);
    const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
    CHECK(result.convergence.residual == compiled.evaluate(result.coordinates).residual_norm());
    const std::size_t components = model.parameters().component_count();
    const double rt = epcsaft::GAS_CONSTANT_J_PER_MOL_K * result.temperature.front();
    const auto charges = epcsaft::equilibrium::species_z(model.parameters());
    for (std::size_t phase = 0; phase < result.phase_count; ++phase) {
        const auto first = result.mole_fractions.begin()
            + static_cast<std::ptrdiff_t>(phase * components);
        const std::vector<double> x(first, first + static_cast<std::ptrdiff_t>(components));
        const double volume = 1.0 / result.molar_densities[phase];
        const auto derivative = model.evaluate(
            result.temperature.front(), volume, x, epcsaft::DerivativeOrder::Gradient);
        double gibbs = derivative.value + result.pressure.front() * volume / rt;
        double lambda_Bx = 0.0;
        double charge = 0.0;
        std::vector<double> chemical(components);
        for (std::size_t species = 0; species < components; ++species) {
            chemical[species] = derivative.gradient[2 + species]
                + std::log(std::max(x[species] / volume,
                                    std::numeric_limits<double>::min()))
                + compiled.standard_gibbs[species] / rt;
            gibbs += x[species] * (
                std::log(std::max(x[species], std::numeric_limits<double>::min()))
                - std::log(volume) - 1.0 + compiled.standard_gibbs[species] / rt);
            charge += charges[species] * x[species];
            for (std::size_t row = 0; row < compiled.rank; ++row)
                lambda_Bx += result.multipliers[row] * compiled.B[row][species] * x[species];
        }
        CHECK(std::abs(gibbs - lambda_Bx) <= limit);
        bool ion_phase = false;
        for (const auto species : compiled.problem.phases[phase].support)
            ion_phase = ion_phase || charges[species] != 0;
        if (ion_phase) CHECK(std::abs(charge) <= limit);
        for (const auto& reaction : compiled.problem.reactions) {
            double affinity = 0.0;
            for (std::size_t species = 0; species < components; ++species)
                affinity += reaction.stoichiometry[species] * chemical[species];
            CHECK(std::abs(affinity) <= limit);
        }
    }
    double feed_dual = 0.0;
    for (std::size_t row = 0; row < compiled.rank; ++row)
        feed_dual += result.multipliers[row] * compiled.balance[row];
    CHECK(std::abs(result.objective - feed_dual) <= limit);
}

void check_s1() {
    epcsaft::HelmholtzModel binary_model(model_support::neutral_binary());
    const auto liquid = equilibrium_support::flash(2);
    epcsaft::core::Capture capture;
    const auto liquid_result = epcsaft::equilibrium::solve_equilibrium(binary_model, liquid);
    verify_s1(binary_model, liquid, liquid_result);
    const auto span = std::find_if(capture.trace().spans.begin(), capture.trace().spans.end(),
        [](const auto& item) { return std::string(item.name) == "solve_equilibrium"; });
    CHECK(span != capture.trace().spans.end() && span->convergence.has_value());
    if (span != capture.trace().spans.end() && span->convergence.has_value()) {
        const auto& record = *span->convergence;
        CHECK(record.residual == liquid_result.convergence.residual);
        CHECK(record.criterion == liquid_result.convergence.criterion);
        CHECK(record.floor == liquid_result.convergence.floor);
        CHECK(record.amplification == liquid_result.convergence.amplification);
        CHECK(record.contraction == liquid_result.convergence.contraction);
        CHECK(record.iterations == liquid_result.convergence.iterations);
        CHECK(record.evaluations == liquid_result.convergence.evaluations);
        CHECK(record.stop == liquid_result.convergence.stop);
    }
    auto vapor = equilibrium_support::flash(2);
    vapor.phases[0].kind = epcsaft::equilibrium::PhaseKind::Vapor;
    verify_s1(binary_model, vapor,
              epcsaft::equilibrium::solve_equilibrium(binary_model, vapor));
    auto flash = equilibrium_support::flash(3);
    flash.phases = {equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor),
                    equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid)};
    flash.temperature.value = 180.0;
    epcsaft::HelmholtzModel ternary_model(model_support::neutral_ternary());
    verify_s1(ternary_model, flash,
              epcsaft::equilibrium::solve_equilibrium(ternary_model, flash));
    epcsaft::HelmholtzModel reactive_model(equilibrium_support::ascani_reactive_ternary());
    const auto reactive = equilibrium_support::ascani_reactive_problem();
    verify_s1(reactive_model, reactive,
              epcsaft::equilibrium::solve_equilibrium(reactive_model, reactive));
    auto ethane = equilibrium_support::flash(1);
    ethane.phases = {equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid),
                     equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor,
                                                epcsaft::equilibrium::AmountSpec::pinned(0.0))};
    ethane.temperature.value = 250.0;
    ethane.pressure = epcsaft::equilibrium::FreeValue::variable(1.0e5);
    epcsaft::HelmholtzModel ethane_model(model_support::ethane());
    verify_s1(ethane_model, ethane,
              epcsaft::equilibrium::solve_equilibrium(ethane_model, ethane));
    auto salt = equilibrium_support::flash(3);
    salt.feed.values = {0.98, 0.01, 0.01};
    salt.temperature.value = 298.15;
    epcsaft::HelmholtzModel salt_model(model_support::aqueous_electrolyte());
    verify_s1(salt_model, salt,
              epcsaft::equilibrium::solve_equilibrium(salt_model, salt));
}

void check_s2() {
    epcsaft::HelmholtzModel model(model_support::neutral_ternary());
    auto binary = equilibrium_support::flash(3);
    binary.phases = {equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor),
                     equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid)};
    binary.temperature.value = 180.0;
    binary.pressure.value = 1.0e5;
    const auto result = epcsaft::equilibrium::solve_equilibrium(model, binary);
    CHECK(result.success && result.reduced_hessian_minimum > 0.0);
    const auto reactive = epcsaft::equilibrium::solve_equilibrium(
        epcsaft::HelmholtzModel(equilibrium_support::ascani_reactive_ternary()),
        equilibrium_support::ascani_reactive_problem());
    CHECK(reactive.success && reactive.reduced_hessian_minimum > 0.0);
}

void check_s6() {
    epcsaft::HelmholtzModel model(equilibrium_support::ascani_reactive_ternary());
    auto two = equilibrium_support::ascani_reactive_problem();
    two.reactions.clear();
    auto one = two;
    one.phases.resize(1);
    const auto one_result = epcsaft::equilibrium::solve_equilibrium(model, one);
    const auto two_result = epcsaft::equilibrium::solve_equilibrium(model, two);
    CHECK(one_result.success && two_result.success);
    if (!one_result.success || !two_result.success) return;
    CHECK(two_result.objective < one_result.objective
          - convergence_limit(two_result.convergence));
}

void check_s7() {
    epcsaft::HelmholtzModel model(equilibrium_support::ascani_reactive_ternary());
    auto problem = equilibrium_support::ascani_reactive_problem();
    problem.reactions.clear();
    const auto first = epcsaft::equilibrium::solve_equilibrium(model, problem);
    CHECK(first.success);
    if (!first.success) return;
    for (const double fraction : {0.3, 0.7}) {
        for (std::size_t i = 0; i < 3; ++i)
            problem.feed.values[i] = (1.0 - fraction) * first.mole_fractions[i]
                + fraction * first.mole_fractions[3 + i];
        const auto result = epcsaft::equilibrium::solve_equilibrium(model, problem, first);
        CHECK(result.success);
        if (!result.success) continue;
        for (std::size_t i = 0; i < 6; ++i)
            CHECK(std::abs(result.mole_fractions[i] - first.mole_fractions[i]) <= 1.0e-8);
    }
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
    epcsaft::HelmholtzModel model(model_support::neutral_ternary());
    auto problem = equilibrium_support::flash(3);
    problem.phases = {equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor),
                      equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid)};
    problem.temperature.value = 180.0;
    const auto solved = epcsaft::equilibrium::solve_equilibrium(model, problem);
    CHECK(solved.success);
    if (!solved.success) return;
    const auto flash = epcsaft::equilibrium::compile_problem(model, problem);
    const auto flash_point = solved.coordinates;
    check_finite_difference(flash, flash_point);
    std::vector<double> primal(flash_point.begin(),
                               flash_point.begin() + static_cast<std::ptrdiff_t>(flash.primal_size));
    epcsaft::equilibrium::GibbsMinimizer nlp(flash, primal);
    Ipopt::Index n = 0, constraints = 0, jacobian_entries = 0, hessian_entries = 0;
    Ipopt::TNLP::IndexStyleEnum style = Ipopt::TNLP::C_STYLE;
    CHECK(nlp.get_nlp_info(n, constraints, jacobian_entries, hessian_entries, style));
    std::vector<Ipopt::Index> rows(static_cast<std::size_t>(hessian_entries));
    std::vector<Ipopt::Index> columns(static_cast<std::size_t>(hessian_entries));
    const std::vector<Ipopt::Number> lambda = solved.constraint_multipliers;
    CHECK(lambda.size() == static_cast<std::size_t>(constraints));
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
        std::vector<Ipopt::Number> jac_plus(static_cast<std::size_t>(constraints * n));
        std::vector<Ipopt::Number> jac_minus(static_cast<std::size_t>(constraints * n));
        CHECK(nlp.eval_jac_g(n, plus.data(), true, constraints, jacobian_entries,
                             nullptr, nullptr, jac_plus.data()));
        CHECK(nlp.eval_jac_g(n, minus.data(), true, constraints, jacobian_entries,
                             nullptr, nullptr, jac_minus.data()));
        for (std::size_t row = column; row < primal.size(); ++row) {
            for (Ipopt::Index constraint = 0; constraint < constraints; ++constraint) {
                grad_plus[row] += lambda[constraint] * jac_plus[constraint * n + row];
                grad_minus[row] += lambda[constraint] * jac_minus[constraint * n + row];
            }
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
    problem.phases.front().composition_guess = {2.0 / 9.0, 4.0 / 9.0, 1.0 / 3.0};
    problem.reactions = reactions;
    double previous = std::numeric_limits<double>::infinity();
    for (const double pressure : {1.0e4, 1.0e3, 1.0e2}) {
        problem.pressure.value = pressure;
        const auto result = epcsaft::equilibrium::solve_equilibrium(model, problem);
        CHECK(result.success && result.mole_fractions.size() == 3);
        if (!result.success || result.mole_fractions.size() != 3) {
            std::fprintf(stderr, "S11 solve: %s\n", result.message.c_str());
            return;
        }
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
    epcsaft::HelmholtzModel model(model_support::neutral_ternary());
    auto problem = equilibrium_support::flash(3);
    problem.phases = {
        equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor),
        equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid),
    };
    problem.temperature.value = 180.0;
    problem.pressure.value = 3.0e6;
    const auto result = epcsaft::equilibrium::solve_equilibrium(model, problem);
    CHECK(!result.success && !result.has_state());
    CHECK(result.message.find("phase") != std::string::npos);
    epcsaft::HelmholtzModel ethane_model(model_support::ethane());
    auto vanished = equilibrium_support::flash(1);
    vanished.phases = {
        equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor,
                                   epcsaft::equilibrium::AmountSpec::pinned(1.0)),
        equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid),
    };
    vanished.temperature.value = 250.0;
    vanished.pressure = epcsaft::equilibrium::FreeValue::variable(1.0e5);
    epcsaft::core::Capture deficit_capture;
    const auto failure = epcsaft::equilibrium::solve_equilibrium(ethane_model, vanished);
    CHECK(!failure.success && !failure.has_state());
    double expected_deficit = 1.0;
    for (std::size_t row = 0; row < failure.rows.size(); ++row)
        if (failure.rows[row].kind == epcsaft::equilibrium::RowKind::Chemical
            && failure.rows[row].phase == 1) expected_deficit -= std::exp(-failure.residuals[row]);
    const auto deficit = std::find_if(deficit_capture.trace().events.begin(),
        deficit_capture.trace().events.end(), [](const auto& event) {
            return std::string(event.name) == "michelsen_deficit" && event.iteration == 1;
        });
    CHECK(deficit != deficit_capture.trace().events.end() && deficit->values.size() == 1);
    if (deficit != deficit_capture.trace().events.end() && deficit->values.size() == 1)
        CHECK(std::abs(deficit->values.front() - expected_deficit) <= 1.0e-12);
    auto symmetric = equilibrium_support::ascani_reactive_problem();
    symmetric.phases[0].composition_guess.clear();
    symmetric.phases[1].composition_guess.clear();
    bool refused = false;
    try { (void)epcsaft::equilibrium::solve_equilibrium(
        epcsaft::HelmholtzModel(equilibrium_support::ascani_reactive_ternary()), symmetric); }
    catch (const std::runtime_error& error) {
        refused = std::string(error.what()).find("start_symmetric") != std::string::npos;
    }
    CHECK(refused);
}

void check_s16() {
    epcsaft::HelmholtzModel model(equilibrium_support::ascani_reactive_ternary());
    auto problem = equilibrium_support::ascani_reactive_problem();
    const auto first = epcsaft::equilibrium::solve_equilibrium(model, problem);
    CHECK(first.success);
    if (!first.success) return;
    for (const double fraction : {0.25, 0.5, 0.75}) {
        for (std::size_t i = 0; i < 3; ++i)
            problem.feed.values[i] = (1.0 - fraction) * first.mole_fractions[i]
                + fraction * first.mole_fractions[3 + i];
        const auto result = epcsaft::equilibrium::solve_equilibrium(model, problem, first);
        CHECK(result.success);
        if (!result.success) continue;
        const double limit = convergence_limit(result.convergence);
        for (std::size_t i = 0; i < 6; ++i)
            CHECK(std::abs(result.mole_fractions[i] - first.mole_fractions[i]) <= limit);
        for (std::size_t phase = 0; phase < 2; ++phase) {
            double affinity = 0.0;
            for (std::size_t row = 0; row < result.rows.size(); ++row)
                if (result.rows[row].kind == epcsaft::equilibrium::RowKind::Chemical
                    && result.rows[row].phase == phase)
                    affinity += problem.reactions[0].stoichiometry[result.rows[row].index]
                        * result.residuals[row];
            CHECK(std::abs(affinity) <= limit);
        }
    }
}

}  // namespace

int main() {
    equilibrium_support::run_check("S1", [] { check_s1(); });
    equilibrium_support::run_check("S2", [] { check_s2(); });
    equilibrium_support::run_check("S6", [] { check_s6(); });
    equilibrium_support::run_check("S7", [] { check_s7(); });
    equilibrium_support::run_check("S8", [] { check_s8(); });
    equilibrium_support::run_check("S10", [] { check_s10(); });
    equilibrium_support::run_check("S11", [] { check_s11(); });
    equilibrium_support::run_check("S13", [] { check_s13(); });
    equilibrium_support::run_check("S15", [] { check_s15(); });
    equilibrium_support::run_check("S16", [] { check_s16(); });
    return equilibrium_support::failures == 0 ? 0 : 1;
}
