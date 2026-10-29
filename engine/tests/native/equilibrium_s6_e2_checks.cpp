// S6-S16, C1, C3, F5 and E2: solve, reuse, availability and thermal-solve checks.
#include "equilibrium_checks_shared.hpp"
#include "one_component_two_density_inputs.hpp"

namespace equilibrium_checks {

void check_s6() {
    epcsaft::HelmholtzModel model(equilibrium_support::synthetic_reactive_ternary());
    auto two = equilibrium_support::synthetic_reactive_problem();
    two.reactions.clear();
    auto one = two;
    one.phases.resize(1);
    const auto one_result = epcsaft::equilibrium::solve_equilibrium(model, one);
    const auto two_result = epcsaft::equilibrium::solve_equilibrium(model, two);
    recording_support::compare(
        "one-success", "success", "bool", one_result.success ? 1.0 : 0.0, 1.0,
        recording_support::fields({recording_support::field("T", 298.15),
                                   recording_support::field("P", 1.0e5),
                                   recording_support::field("feed", two.feed.values),
                                   recording_support::field("phase_count", 1U)}),
        "Relationship", "solve", "one-phase solve succeeds", "reference_or_one", 0.0, 0.0,
        one_result.success, "s6-one-success", "one", "success", "T,P,feed,phases", "success=1", "case");
    recording_support::compare(
        "two-success", "success", "bool", two_result.success ? 1.0 : 0.0, 1.0,
        recording_support::fields({recording_support::field("T", 298.15),
                                   recording_support::field("P", 1.0e5),
                                   recording_support::field("feed", two.feed.values),
                                   recording_support::field("phase_count", 2U)}),
        "Relationship", "solve", "two-phase solve succeeds", "reference_or_one", 0.0, 0.0,
        two_result.success, "s6-two-success", "two", "success", "T,P,feed,phases", "success=1", "case");
    const double objective_gap = one_result.objective - two_result.objective;
    const double objective_limit = two_result.convergence.criterion == nullptr
        ? 0.0 : convergence_limit(two_result.convergence);
    const bool objective_gap_ok = std::isfinite(objective_gap) && std::isfinite(objective_limit)
        && objective_gap > objective_limit;
    recording_support::compare_custom(
        "objective-gap", "objective_gap", "1", objective_gap, objective_limit,
        recording_support::fields({recording_support::field("T", 298.15),
                                   recording_support::field("P", 1.0e5),
                                   recording_support::field("feed", two.feed.values),
                                   recording_support::field("one_objective", one_result.objective),
                                   recording_support::field("two_objective", two_result.objective),
                                   recording_support::field("convergence_limit", objective_limit)}),
        "Relationship", "phase-split", "two-phase objective gap exceeds convergence limit",
        "isfinite(model) && isfinite(reference) && model > reference", objective_gap_ok,
        "s6-objective-gap", "two-minus-one", "\\Delta g", "T,P,feed,objective,limit",
        "\\Delta g>\\ell", "objective");
    CHECK(one_result.success && two_result.success);
    if (!one_result.success || !two_result.success) return;
    CHECK(two_result.objective < one_result.objective
          - convergence_limit(two_result.convergence));
}

void check_s7() {
    epcsaft::HelmholtzModel model(equilibrium_support::synthetic_reactive_ternary());
    auto problem = equilibrium_support::synthetic_reactive_problem();
    problem.reactions.clear();
    const auto first = epcsaft::equilibrium::solve_equilibrium(model, problem);
    CHECK(first.success);
    for (const double fraction : {0.3, 0.7}) {
        for (std::size_t i = 0; i < 3; ++i)
            problem.feed.values[i] = (1.0 - fraction) * first.mole_fractions[i]
                + fraction * first.mole_fractions[3 + i];
        const auto result = epcsaft::equilibrium::solve_equilibrium(model, problem, first);
        CHECK(result.success);
        if (!result.success) continue;
        for (std::size_t i = 0; i < 6; ++i) {
            const double difference = std::abs(
                result.mole_fractions[i] - first.mole_fractions[i]);
            const auto point = recording_support::number(fraction) + "/" + std::to_string(i);
            recording_support::compare(
                (std::string("continuation/") + point).c_str(),
                "mole_fraction_difference", "1", difference, 0.0,
                recording_support::fields({recording_support::field("T", problem.temperature.value),
                                           recording_support::field("P", problem.pressure.value),
                                           recording_support::field("feed", problem.feed.values),
                                           recording_support::field("fraction", fraction),
                                           recording_support::field("component", i),
                                           recording_support::field("reference_solution", first.mole_fractions),
                                           recording_support::field("solution", result.mole_fractions)}),
                "Relationship", "continuation", "continuation preserves the solved composition",
                "reference_or_one", 1.0e-8, 0.0, difference <= 1.0e-8,
                "s7-continuation", point.c_str(), "\\Delta x_i",
                "T,P,feed,\\alpha,component", "\\Delta x_i\\leq10^{-8}", "component");
        }
    }
}

void check_c1_multiplier_scale() {
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    auto small = equilibrium_support::flash(2);
    auto large = small;
    large.feed.values = {5.0, 5.0};
    const auto small_result = epcsaft::equilibrium::solve_equilibrium(model, small);
    const auto large_result = epcsaft::equilibrium::solve_equilibrium(model, large);
    CHECK(small_result.success && large_result.success);
    if (!small_result.success || !large_result.success) return;

    const bool same_intensives = small_result.temperature == large_result.temperature
        && small_result.pressure == large_result.pressure;
    const bool same_amount_shape = small_result.amounts.size() == large_result.amounts.size();
    if (!same_amount_shape) {
        CHECK(same_intensives);
        CHECK(same_amount_shape);
        return;
    }
    for (std::size_t i = 0; i < small_result.amounts.size(); ++i) {
        const double ratio = large_result.amounts[i] / small_result.amounts[i];
        recording_support::compare(
            (std::string("amount-scale/") + std::to_string(i)).c_str(), "phase_amount_ratio", "1",
            ratio, 10.0,
            recording_support::fields({recording_support::field("T", *small_result.temperature),
                                       recording_support::field("P", *small_result.pressure),
                                       recording_support::field("feed_small", small.feed.values),
                                       recording_support::field("feed_large", large.feed.values),
                                       recording_support::field("phase", i)}),
            "Relationship", "extensive-scaling", "amounts scale with feed by ten",
            "strict_absolute", 1.0e-7, 0.0, std::abs(ratio - 10.0) < 1.0e-7,
        "c1-multiplier-scale", std::to_string(i).c_str(), "n", "T,P,feed,phase",
        "n_{large}/n_{small}=10", "phase");
    }
    CHECK(same_intensives);
    CHECK(same_amount_shape);

    const auto small_compiled = epcsaft::equilibrium::compile_problem(model, small);
    const auto large_compiled = epcsaft::equilibrium::compile_problem(model, large);
    const auto small_rows = callback_rows(small_compiled);
    const auto large_rows = callback_rows(large_compiled);
    CHECK(small_rows == large_rows);
    CHECK(small_result.constraint_multipliers.size() == small_rows.size());
    CHECK(large_result.constraint_multipliers.size() == large_rows.size());
    if (small_rows != large_rows
        || small_result.constraint_multipliers.size() != small_rows.size()
        || large_result.constraint_multipliers.size() != large_rows.size()) {
        return;
    }

    const std::vector<double> small_point(
        small_result.coordinates.begin(),
        small_result.coordinates.begin()
            + static_cast<std::ptrdiff_t>(small_compiled.primal_size));
    const std::vector<double> large_point(
        large_result.coordinates.begin(),
        large_result.coordinates.begin()
            + static_cast<std::ptrdiff_t>(large_compiled.primal_size));
    const auto small_objective = small_compiled.objective_gradient(small_point);
    const auto large_objective = large_compiled.objective_gradient(large_point);
    const auto small_kkt = weighted_gradient(
        small_compiled, small_point, small_result.constraint_multipliers, 1.0);
    const auto large_kkt = weighted_gradient(
        large_compiled, large_point, large_result.constraint_multipliers, 1.0);
    const double small_kkt_scale = std::max(1.0, maximum_abs(small_objective));
    const double large_kkt_scale = std::max(1.0, maximum_abs(large_objective));
    const double small_kkt_limit =
        100.0 * epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value * small_kkt_scale;
    const double large_kkt_limit =
        100.0 * epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value * large_kkt_scale;
    const double small_kkt_residual = maximum_abs(small_kkt);
    const double large_kkt_residual = maximum_abs(large_kkt);
    const auto kkt_inputs = recording_support::fields({
        recording_support::field("T", 300.0),
        recording_support::field("P", 1.0e5),
        recording_support::field("feed_small", small.feed.values),
        recording_support::field("feed_large", large.feed.values),
        recording_support::field("feed_total_small", small_compiled.feed_total),
        recording_support::field("feed_total_large", large_compiled.feed_total),
        recording_support::field("small_solver_status", small_result.solver_status),
        recording_support::field("large_solver_status", large_result.solver_status),
        recording_support::field("small_convergence_stop",
            epcsaft::core::stop_name(small_result.convergence.stop)),
        recording_support::field("large_convergence_stop",
            epcsaft::core::stop_name(large_result.convergence.stop)),
        recording_support::field("active_constraint_count", small_rows.size()),
        recording_support::field("small_objective_gradient", small_objective),
        recording_support::field("large_objective_gradient", large_objective),
        recording_support::field("small_constraint_multipliers",
            small_result.constraint_multipliers),
        recording_support::field("large_constraint_multipliers",
            large_result.constraint_multipliers),
        recording_support::field("small_kkt_scale", small_kkt_scale),
        recording_support::field("large_kkt_scale", large_kkt_scale),
    });
    // Since the #87 correction the returned constraint multipliers are formed at the returned point
    // from the row multipliers, so these rows verify that reconstruction's KKT residual and feed
    // scaling, no longer the optimizer's own duals.
    recording_support::compare(
        "c1-multiplier-scale/kkt-small", "kkt_residual", "1", small_kkt_residual, 0.0,
        kkt_inputs, "Relationship", "KKT", "small-feed returned multipliers satisfy KKT",
        "reference_or_one", small_kkt_limit, 0.0,
        small_kkt_residual <= small_kkt_limit, "c1-multiplier-scale", "small",
        "||r_KKT||", "T,P,feed,lambda,scale", "||r_KKT||<=100*row_residual*scale", "feed");
    recording_support::compare(
        "c1-multiplier-scale/kkt-large", "kkt_residual", "1", large_kkt_residual, 0.0,
        kkt_inputs, "Relationship", "KKT", "large-feed returned multipliers satisfy KKT",
        "reference_or_one", large_kkt_limit, 0.0,
        large_kkt_residual <= large_kkt_limit, "c1-multiplier-scale", "large",
        "||r_KKT||", "T,P,feed,lambda,scale", "||r_KKT||<=100*row_residual*scale", "feed");
    CHECK(small_kkt_residual <= small_kkt_limit);
    CHECK(large_kkt_residual <= large_kkt_limit);

    const double feed_scale = large_compiled.feed_total / small_compiled.feed_total;
    double multiplier_difference = 0.0;
    for (std::size_t i = 0; i < small_result.constraint_multipliers.size(); ++i) {
        const double expected = feed_scale * small_result.constraint_multipliers[i];
        multiplier_difference = std::max(
            multiplier_difference, std::abs(large_result.constraint_multipliers[i] - expected)
                / std::max(1.0, std::abs(expected)));
    }
    recording_support::compare(
        "c1-multiplier-scale/multipliers", "relative_multiplier_difference", "1",
        multiplier_difference, 0.0, kkt_inputs, "Relationship", "extensive-scaling",
        "constraint multipliers scale with feed amount", "reference_or_one", 1.0e-7, 0.0,
        multiplier_difference <= 1.0e-7, "c1-multiplier-scale", "multipliers",
        "delta_lambda", "T,P,feed,lambda", "lambda_large=10*lambda_small", "feed");

    double intensive_difference = 0.0;
    intensive_difference = std::max(
        intensive_difference,
        std::abs(*small_result.temperature - *large_result.temperature)
            / std::max(1.0, std::abs(*small_result.temperature)));
    intensive_difference = std::max(
        intensive_difference,
        std::abs(*small_result.pressure - *large_result.pressure)
            / std::max(1.0, std::abs(*small_result.pressure)));
    for (std::size_t i = 0; i < small_result.mole_fractions.size(); ++i) {
        intensive_difference = std::max(
            intensive_difference,
            std::abs(small_result.mole_fractions[i] - large_result.mole_fractions[i]));
    }
    for (std::size_t i = 0; i < small_result.molar_densities.size(); ++i) {
        intensive_difference = std::max(
            intensive_difference,
            std::abs(small_result.molar_densities[i] - large_result.molar_densities[i])
                / std::max(1.0, std::abs(small_result.molar_densities[i])));
        intensive_difference = std::max(
            intensive_difference,
            std::abs(small_result.packing_fractions[i] - large_result.packing_fractions[i]));
    }
    recording_support::compare(
        "c1-multiplier-scale/intensives", "relative_intensive_difference", "1",
        intensive_difference, 0.0, kkt_inputs, "Relationship", "extensive-scaling",
        "intensive state is invariant under feed scaling", "reference_or_one", 1.0e-8, 0.0,
        intensive_difference <= 1.0e-8, "c1-multiplier-scale", "intensives",
        "delta_state", "T,P,feed", "state_large=state_small", "feed");
    CHECK(small_result.reduced_hessian_minimum > 0.0);
    CHECK(large_result.reduced_hessian_minimum > 0.0);

    auto double_scaled = large_result.constraint_multipliers;
    for (double& value : double_scaled) value *= feed_scale;
    const double double_scaled_kkt = maximum_abs(
        weighted_gradient(large_compiled, large_point, double_scaled, 1.0));
    recording_support::compare_custom(
        "c1-multiplier-scale/double-scaling", "kkt_residual", "1", double_scaled_kkt,
        large_kkt_limit, kkt_inputs, "Relationship", "KKT-regression-control",
        "extra feed scaling violates the KKT residual", "model > reference",
        double_scaled_kkt > large_kkt_limit, "c1-multiplier-scale", "double-scaling",
        "||r_KKT||", "T,P,feed,lambda", "double-scaled lambda rejected", "feed");

    auto tiny = small;
    tiny.feed.values = {5.0e-10, 5.0e-10};
    const auto tiny_result = epcsaft::equilibrium::solve_equilibrium(model, tiny);
    CHECK(tiny_result.success);
    const auto tiny_compiled = epcsaft::equilibrium::compile_problem(model, tiny);
    const auto tiny_rows = callback_rows(tiny_compiled);
    CHECK(tiny_rows == small_rows);
    CHECK(tiny_result.constraint_multipliers.size() == tiny_rows.size());
    if (tiny_rows != small_rows
        || tiny_result.constraint_multipliers.size() != tiny_rows.size()) return;

    const auto tiny_point = std::vector<double>(
        tiny_result.coordinates.begin(),
        tiny_result.coordinates.begin()
            + static_cast<std::ptrdiff_t>(tiny_compiled.primal_size));
    const auto tiny_objective = tiny_compiled.objective_gradient(tiny_point);
    const auto tiny_kkt = weighted_gradient(
        tiny_compiled, tiny_point, tiny_result.constraint_multipliers, 1.0);
    const double tiny_kkt_scale = maximum_abs(tiny_objective);
    CHECK(tiny_kkt_scale > 0.0);
    const double tiny_kkt_limit =
        100.0 * epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value * tiny_kkt_scale;
    const double tiny_kkt_residual = maximum_abs(tiny_kkt);
    const auto tiny_inputs = recording_support::fields({
        recording_support::field("T", 300.0),
        recording_support::field("P", 1.0e5),
        recording_support::field("feed", tiny.feed.values),
        recording_support::field("feed_total", tiny_compiled.feed_total),
        recording_support::field("solver_status", tiny_result.solver_status),
        recording_support::field("convergence_stop",
            epcsaft::core::stop_name(tiny_result.convergence.stop)),
        recording_support::field("residual", tiny_result.convergence.residual),
        recording_support::field("convergence_limit", convergence_limit(tiny_result.convergence)),
        recording_support::field("objective_gradient", tiny_objective),
        recording_support::field("constraint_multipliers",
            tiny_result.constraint_multipliers),
        recording_support::field("kkt_scale", tiny_kkt_scale),
        recording_support::field("kkt_limit", tiny_kkt_limit),
        recording_support::field("reduced_hessian_minimum",
            tiny_result.reduced_hessian_minimum),
    });
    const bool tiny_stop_ok = tiny_result.convergence.stop == epcsaft::core::Stop::converged;
    recording_support::compare(
        "c1-multiplier-scale/tiny-convergence", "convergence_stop", "bool",
        tiny_stop_ok ? 1.0 : 0.0, 1.0, tiny_inputs, "Relationship", "equilibrium-residual",
        "tiny-step is accepted only after residual convergence", "absolute", 0.0, 0.0,
        tiny_stop_ok, "c1-multiplier-scale", "tiny-convergence", "stop", "feed,residual",
        "stop=converged", "feed");
    recording_support::compare(
        "c1-multiplier-scale/tiny-kkt", "kkt_residual", "1", tiny_kkt_residual, 0.0,
        tiny_inputs, "Relationship", "KKT", "tiny-feed returned multipliers satisfy KKT",
        "strict_absolute", tiny_kkt_limit, 0.0, tiny_kkt_residual <= tiny_kkt_limit,
        "c1-multiplier-scale", "tiny", "||r_KKT||", "T,P,feed,lambda,scale",
        "||r_KKT||<=100*row_residual*||grad f||", "feed");
    const bool tiny_reduced_hessian_ok =
        std::isfinite(tiny_result.reduced_hessian_minimum)
        && tiny_result.reduced_hessian_minimum > 0.0;
    recording_support::compare_custom(
        "c1-multiplier-scale/tiny-reduced-hessian", "reduced_hessian_minimum", "1",
        tiny_result.reduced_hessian_minimum, 0.0, tiny_inputs, "Relationship", "KKT",
        "tiny-feed reduced Hessian remains positive", "model > 0", tiny_reduced_hessian_ok,
        "c1-multiplier-scale", "tiny-reduced-hessian", "lambda_min", "feed,lambda",
        "lambda_min>0", "feed");
    CHECK(tiny_stop_ok);
    CHECK(tiny_kkt_residual <= tiny_kkt_limit);
    CHECK(tiny_reduced_hessian_ok);

    const double tiny_feed_scale = tiny_compiled.feed_total / small_compiled.feed_total;
    double tiny_multiplier_difference = 0.0;
    for (std::size_t i = 0; i < tiny_result.constraint_multipliers.size(); ++i) {
        const double expected = tiny_feed_scale * small_result.constraint_multipliers[i];
        const double denominator = std::max(std::abs(expected), tiny_kkt_limit);
        tiny_multiplier_difference = std::max(
            tiny_multiplier_difference,
            std::abs(tiny_result.constraint_multipliers[i] - expected) / denominator);
    }
    recording_support::compare(
        "c1-multiplier-scale/tiny-multipliers", "relative_multiplier_difference", "1",
        tiny_multiplier_difference, 0.0, tiny_inputs, "Relationship", "extensive-scaling",
        "tiny-feed multipliers scale with feed amount", "reference_or_one", 1.0e-6, 0.0,
        tiny_multiplier_difference <= 1.0e-6, "c1-multiplier-scale", "tiny-multipliers",
        "delta_lambda", "feed,lambda,scale", "lambda_tiny=1e-9*lambda_small", "feed");

    double tiny_intensive_difference = 0.0;
    tiny_intensive_difference = std::max(
        tiny_intensive_difference,
        std::abs(*tiny_result.temperature - *small_result.temperature)
            / std::max(1.0, std::abs(*small_result.temperature)));
    tiny_intensive_difference = std::max(
        tiny_intensive_difference,
        std::abs(*tiny_result.pressure - *small_result.pressure)
            / std::max(1.0, std::abs(*small_result.pressure)));
    for (std::size_t i = 0; i < small_result.mole_fractions.size(); ++i) {
        tiny_intensive_difference = std::max(
            tiny_intensive_difference,
            std::abs(tiny_result.mole_fractions[i] - small_result.mole_fractions[i]));
    }
    for (std::size_t i = 0; i < small_result.molar_densities.size(); ++i) {
        tiny_intensive_difference = std::max(
            tiny_intensive_difference,
            std::abs(tiny_result.molar_densities[i] - small_result.molar_densities[i])
                / std::max(1.0, std::abs(small_result.molar_densities[i])));
        tiny_intensive_difference = std::max(
            tiny_intensive_difference,
            std::abs(tiny_result.packing_fractions[i] - small_result.packing_fractions[i]));
    }
    recording_support::compare(
        "c1-multiplier-scale/tiny-intensives", "relative_intensive_difference", "1",
        tiny_intensive_difference, 0.0, tiny_inputs, "Relationship", "extensive-scaling",
        "tiny-feed intensive state is invariant", "reference_or_one", 1.0e-8, 0.0,
        tiny_intensive_difference <= 1.0e-8, "c1-multiplier-scale", "tiny-intensives",
        "delta_state", "feed,T,P", "state_tiny=state_small", "feed");
}

void check_c1_ideal_caller_chart() {
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    auto problem = equilibrium_support::flash(2);
    problem.feed.values = {0.2, 0.9};
    const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
    auto point = compiled.generic_point();
    CHECK(point.size() >= 2);
    if (point.size() < 2) return;
    point[0] = std::log(0.2);
    point[1] = std::log(0.9);
    const auto base = compiled.evaluate(point);
    const auto analytic = compiled.objective_gradient(point);
    const auto inputs = [&](std::size_t column, double step) {
        return recording_support::fields({
            recording_support::field("T", problem.temperature.value),
            recording_support::field("P", problem.pressure.value),
            recording_support::field("feed", problem.feed.values),
            recording_support::field("x", std::vector<double>{0.2, 0.9}),
            recording_support::field("sum_x", 1.1),
            recording_support::field("coordinate", column),
            recording_support::field("step", step),
            recording_support::field("objective", base.objective),
        });
    };
    for (const std::size_t column : {std::size_t(0), std::size_t(1)}) {
        for (const double step : {1.0e-2, 3.0e-3, 1.0e-3}) {
            auto plus = point;
            auto minus = point;
            plus[column] += step;
            minus[column] -= step;
            const double finite_difference =
                (compiled.evaluate(plus).objective - compiled.evaluate(minus).objective)
                / (2.0 * step);
            const double tolerance = step == 1.0e-2 ? 1.0e-4
                : step == 3.0e-3 ? 1.0e-5 : 1.0e-6;
            const bool verdict = close_with_tolerance(analytic[column], finite_difference, tolerance);
            const std::string point_id = std::to_string(column) + "/" + recording_support::number(step);
            recording_support::compare(
                (std::string("objective-gradient/") + point_id).c_str(), "objective_gradient", "1",
                analytic[column], finite_difference, inputs(column, step), "Relationship",
                "centered_objective_value_ladder", "off-normalized caller gradient matches objective",
                "reference_or_one", 0.0, tolerance, verdict, "c1-equilibrium-caller", point_id.c_str(),
                "dG/dz", "T,P,feed,x,sum_x,coordinate,step", "dG/dz=\u0394G/(2\u0394z)", "coordinate,step");
        }
    }
}

void check_f5() {
    epcsaft::HelmholtzModel model(properties_support::neutral_ternary());
    auto problem = equilibrium_support::flash(3);
    problem.phases = {equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor),
                      equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid)};
    problem.temperature.value = 180.0;
    const auto solved = epcsaft::equilibrium::solve_equilibrium(model, problem);
    CHECK(solved.success);
    const auto flash = epcsaft::equilibrium::compile_problem(model, problem);
    const auto flash_point = solved.coordinates;
    check_finite_difference(flash, flash_point);
    std::vector<double> primal(flash_point.begin(),
                               flash_point.begin() + static_cast<std::ptrdiff_t>(flash.primal_size));
    const auto active_rows = callback_rows(flash);
    std::vector<double> lambda = solved.constraint_multipliers;
    CHECK(lambda.size() == active_rows.size());
    if (lambda.size() != active_rows.size()) return;
    for (std::size_t row = 0; row < lambda.size(); ++row) {
        if (std::abs(lambda[row]) < 1.0e-12)
            lambda[row] = 0.25 + 0.1 * static_cast<double>(row);
    }
    // F5 previous-failure check: vary objective factors and nonzero callback
    // multipliers independently, then differentiate the complete weighted gradient.
    auto pinned_problem = equilibrium_support::boundary();
    pinned_problem.temperature.value = 180.0;
    epcsaft::HelmholtzModel pinned_model(properties_support::neutral_binary());
    const auto pinned = epcsaft::equilibrium::compile_problem(
        pinned_model, pinned_problem
    );
    const auto pinned_point = pinned.generic_point();
    const auto pinned_rows = callback_rows(pinned);
    std::vector<double> pinned_lambda(pinned_rows.size());
    for (std::size_t row = 0; row < pinned_lambda.size(); ++row)
        pinned_lambda[row] = 0.35 + 0.1 * static_cast<double>(row);
    check_weighted_hessian_gradient_ladder(flash, primal, lambda, "regular", 1.0e-6);
    check_weighted_hessian_gradient_ladder(pinned, pinned_point, pinned_lambda, "pinned", 1.0e-5);
    // Free intensive values must be differentiated at the callback point,
    // including reaction offsets and simultaneous temperature/pressure freedom.
    for (int variant = 0; variant < 6; ++variant) {
        auto request = pinned_problem;
        request.temperature = epcsaft::equilibrium::FreeValue::variable(180.0);
        request.pressure = epcsaft::equilibrium::FreeValue::fixed(100000.0);
        const char* label = "pinned-free-T";
        if (variant == 1) {
            request.pressure = epcsaft::equilibrium::FreeValue::variable(100000.0);
            request.phases[1].amount = epcsaft::equilibrium::AmountSpec::pinned(0.2);
            auto third_phase = request.phases[1];
            third_phase.name = "L2";
            third_phase.kind = epcsaft::equilibrium::PhaseKind::Liquid;
            third_phase.amount = epcsaft::equilibrium::AmountSpec::pinned(0.3);
            request.phases.push_back(third_phase);
            label = "pinned-free-TP";
        } else if (variant == 2) {
            request.reactions = {{"A-to-B", {-1.0, 1.0}, std::nullopt, 850.0}};
            label = "pinned-reactive-T";
        } else if (variant == 5) {
            // Non-affine correlated standard state: its exact s', s'' enter the T column and
            // curvature; the callback point below evaluates at 313 K inside [250, 450] K.
            request.temperature.value = 300.0;
            request.reactions = {{"A-to-B", {-1.0, 1.0}, std::nullopt, std::nullopt, legacy_isomerisation(1.0, moderate_isomerisation_a())}};
            label = "pinned-correlated-T";
        } else if (variant >= 3) {
            request.temperature.value = 300.0;
            label = variant == 3 ? "pinned-association-T" : "pinned-ionic-T";
            if (variant == 4) request.feed.values = {0.98, 0.01, 0.01};
        }
        epcsaft::HelmholtzModel action_model(variant == 3
            ? properties_support::associating_binary() : variant == 4
            ? properties_support::aqueous_electrolyte() : properties_support::neutral_binary());
        const auto compiled = epcsaft::equilibrium::compile_problem(action_model, request);
        auto callback_point = compiled.generic_point();
        callback_point[compiled.primal_size - compiled.free_intensive_count] = variant >= 3 ? 313.0 : 193.0;
        if (variant == 1) callback_point[compiled.primal_size - 1] = 120000.0;
        std::vector<double> weights(callback_rows(compiled).size(), 0.0);
        for (std::size_t row = 0; row < weights.size(); ++row)
            weights[row] = 0.35 + 0.1 * static_cast<double>(row);
        check_weighted_hessian_gradient_ladder(compiled, callback_point, weights, label, 1.0e-5);
        if (variant == 4) {
            // #641: a charged molecule's neutral d(T) reaches the T column through packing_terms' hand-written dd/dT.
            auto charged = properties_support::aqueous_electrolyte();
            charged.molecular_ion = {0, 1, 0};
            const epcsaft::HelmholtzModel charged_model(charged);
            check_finite_difference(epcsaft::equilibrium::compile_problem(charged_model, request), callback_point);
        }
        if (variant == 1) {
            epcsaft::equilibrium::GibbsMinimizer nlp(compiled, callback_point);
            Ipopt::Index variables = 0, constraints = 0, jacobian_entries = 0, hessian_entries = 0;
            Ipopt::TNLP::IndexStyleEnum index_style;
            CHECK(nlp.get_nlp_info(variables, constraints, jacobian_entries, hessian_entries, index_style));
            const auto rows = callback_rows(compiled);
            std::vector<double> pressure_weights(rows.size(), 0.0);
            const auto phase_point = compiled.prepare(callback_point);
            std::size_t selected = 0;
            while (compiled.rows[rows[selected]].kind != epcsaft::equilibrium::RowKind::Pressure) ++selected;
            pressure_weights[selected] = 0.7;
            const std::size_t phase = compiled.rows[rows[selected]].phase;
            const std::size_t w = phase_point.shape[phase].volume;
            const std::size_t P = compiled.primal_size - 1;
            std::vector<double> hessian(static_cast<std::size_t>(hessian_entries));
            CHECK(nlp.eval_h(variables, callback_point.data(), true, 0.0, constraints,
                pressure_weights.data(), true, hessian_entries, nullptr, nullptr, hessian.data()));
            const double reference = pressure_weights[selected] * phase_point.volume[phase]
                / (epcsaft::GAS_CONSTANT_J_PER_MOL_K * phase_point.temperature) * compiled.problem.pressure.value;
            const double action = hessian[P * (P + 1) / 2 + w] * compiled.problem.pressure.value;
            const double allowance = 1.0e-12 + 1.0e-8 * std::abs(reference);
            const bool agrees = std::abs(action - reference) <= allowance && std::abs(reference) > allowance;
            recording_support::compare("hessian/pressure-packing-identity", "P_scale_H_Pw", "1",
                action, reference, recording_support::fields({recording_support::field("V_m3_per_mol", phase_point.volume[phase]),
                    recording_support::field("T_K", phase_point.temperature), recording_support::field("lambda", 0.7),
                    recording_support::field("P_scale_Pa", compiled.problem.pressure.value)}),
                "Relationship", "pressure-row-identity", "H_Pw = lambda V/(RT)",
                "reference_scale", 1.0e-12, 1.0e-8, agrees);
        }
        if (variant == 5) {
            // The pinned ladder criterion (1e-5) is looser than the affine-shortcut defect at this
            // point (~2e-6), so the T,T curvature entry is compared directly at 1e-9 scaled: the
            // falsifier flips if the retired s'' = -2 s'/T shortcut returns.
            epcsaft::equilibrium::GibbsMinimizer nlp(compiled, callback_point);
            Ipopt::Index variables = 0, constraints = 0, jacobian_entries = 0, hessian_entries = 0;
            Ipopt::TNLP::IndexStyleEnum index_style;
            CHECK(nlp.get_nlp_info(variables, constraints, jacobian_entries, hessian_entries, index_style));
            const std::size_t t_index = compiled.primal_size - compiled.free_intensive_count;
            std::vector<Ipopt::Number> hessian(static_cast<std::size_t>(hessian_entries), 0.0);
            CHECK(nlp.eval_h(variables, callback_point.data(), true, 0.0, constraints,
                             weights.data(), true, hessian_entries, nullptr, nullptr, hessian.data()));
            const double h_tt = hessian[t_index * (t_index + 1) / 2 + t_index];
            const double tt_step = centered_step(callback_point[t_index]);
            auto plus = callback_point, minus = callback_point;
            plus[t_index] += tt_step; minus[t_index] -= tt_step;
            const double h_tt_centered = (weighted_gradient(compiled, plus, weights, 0.0)[t_index]
                - weighted_gradient(compiled, minus, weights, 0.0)[t_index]) / (2.0 * tt_step);
            const bool tt_ok = close_with_tolerance(h_tt, h_tt_centered, 1.0e-9);
            recording_support::compare("correlated-jacobian/hessian-TT", "weighted_lagrangian_hessian", "coordinate",
                h_tt, h_tt_centered, recording_support::fields({equilibrium_support::weighted_hessian_context(compiled),
                    recording_support::field("point", callback_point), recording_support::field("step", tt_step)}),
                "Relationship", "weighted_hessian_gradient_finite_difference",
                "correlated standard-state T,T curvature is exact beyond the affine shortcut",
                "reference_or_one", 0.0, 1.0e-9, tt_ok, "f5-correlated-jacobian", "hessian-TT",
                "\\nabla^2\\mathcal L_{TT}", "\\mathbf{z},\\lambda,\\Delta T", "exact s''", "entry");
            // Rows, Jacobian temperature column and lever column against centered differences.
            std::vector<double> full = callback_point;
            full.resize(compiled.unknown_count(), 0.0);
            const auto jacobian = compiled.jacobian(full);
            const auto levers = compiled.lever_jacobian(full);
            const std::size_t lever_width = levers.size() / compiled.rows.size();
            const std::size_t column = compiled.primal_size - compiled.free_intensive_count;
            const double step = 1.0e-3;
            auto warmer = full, colder = full;
            warmer[column] += step; colder[column] -= step;
            const auto warm = compiled.evaluate(warmer).values, cold = compiled.evaluate(colder).values;
            for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
                const double centered = (warm[row] - cold[row]) / (2.0 * step);
                const double analytic = jacobian[row * full.size() + column];
                const bool verdict = close(analytic, centered)
                    && (compiled.rows[row].kind != epcsaft::equilibrium::RowKind::Chemical
                        || close(levers[row * lever_width], centered));
                const auto id = std::string("correlated-jacobian/") + std::to_string(row);
                recording_support::compare(id.c_str(), "equilibrium_jacobian", "1", analytic, centered,
                    recording_support::fields({equilibrium_support::weighted_hessian_context(compiled),
                        recording_support::field("point", full), recording_support::field("row", row),
                        recording_support::field("lever_column", levers[row * lever_width]),
                        recording_support::field("delta_T", step)}),
                    "Relationship", "finite-difference",
                    "correlated standard-state temperature column matches centered rows",
                    "reference_or_one", 0.0, 1.0e-6, verdict, "f5-correlated-jacobian", id.c_str(),
                    "J_{r,T}", "T,P,\\mathbf{x},\\Delta T", "J_{r,T}=\\Delta r/\\Delta T", "row");
            }
        }
    }
}

void check_s11() {
    const std::vector<epcsaft::equilibrium::ReactionSpec> reactions = {
        {"isomerisation", {-1.0, 1.0, 0.0}, 2.0, std::nullopt}
    };
    const auto standard = epcsaft::equilibrium::standard_state(300.0, 3, reactions, {});
    const auto affinity = equilibrium_support::reaction_affinities(reactions, standard);
    const double affinity_reference = -epcsaft::GAS_CONSTANT_J_PER_MOL_K * 300.0 * std::log(2.0);
    recording_support::compare(
        "standard-affinity", "reaction_affinity", "J/mol", affinity.front(), affinity_reference,
        recording_support::fields({recording_support::field("T", 300.0),
                                   recording_support::field("stoichiometry", std::vector<double>{-1.0, 1.0, 0.0}),
                                   recording_support::field("K", 2.0)}),
        "Relationship", "reaction-standard-state", "standard affinity follows K=2",
        "strict_absolute", 1e-10, 0.0,
        std::abs(affinity.front() - affinity_reference) < 1e-10,
        "s11-standard-affinity", "reaction", "A", "T,\\nu,K", "A=-RT\\ln K", "reaction");
    CHECK(affinity.size() == 1);
    CHECK(std::isfinite(affinity.front()));
    CHECK(std::abs(affinity.front() + epcsaft::GAS_CONSTANT_J_PER_MOL_K * 300.0
                   * std::log(2.0)) < 1.0e-10);
    epcsaft::HelmholtzModel model(equilibrium_support::synthetic_reactive_ternary());
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
            + affinity.front() / (epcsaft::GAS_CONSTANT_J_PER_MOL_K * *result.temperature));
        const auto pressure_id = recording_support::number(pressure);
        const auto previous_gap = std::isfinite(previous)
            ? recording_support::number(previous) : std::string("none");
        const auto gap_inputs = recording_support::fields({
            recording_support::field("T", *result.temperature),
            recording_support::field("P", pressure),
            recording_support::field("x", result.mole_fractions),
            recording_support::field("affinity", affinity.front()),
            recording_support::field("previous_gap", previous_gap),
        });
        const bool gap_contracts = std::isfinite(gap) && std::isfinite(previous)
            && gap <= previous;
        if (std::isfinite(previous)) {
            recording_support::compare_custom(
                (std::string("pressure-gap/") + pressure_id).c_str(),
                "reaction_equilibrium_gap", "1", gap, previous, gap_inputs,
                "Relationship", "reaction-equilibrium",
                "reaction gap contracts as pressure decreases",
                "isfinite(model) && isfinite(reference) && model <= reference", gap_contracts,
                "s11-gap", pressure_id.c_str(), "\\Delta", "T,P,\\mathbf{x},A,\\Delta_{previous}",
                "\\Delta(P)\\leq\\Delta(P_{previous})", "pressure");
        } else {
            recording_support::observe(
                (std::string("pressure-gap/") + pressure_id).c_str(),
                "reaction_equilibrium_gap", "1", gap, std::nullopt, gap_inputs,
                "Relationship", "reaction-equilibrium",
                "first reaction gap establishes the contraction baseline",
                "s11-gap", pressure_id.c_str(), "\\Delta", "T,P,\\mathbf{x},A",
                "\\Delta(P_{first})", "pressure");
        }
        CHECK(std::isfinite(gap) && (std::isinf(previous) || gap_contracts));
        previous = gap;
        if (pressure == 1.0e2) {
            recording_support::compare(
                "final-gap", "reaction_equilibrium_gap", "1", gap, 0.0,
                recording_support::fields({recording_support::field("T", *result.temperature),
                                           recording_support::field("P", pressure),
                                           recording_support::field("x", result.mole_fractions),
                                           recording_support::field("affinity", affinity.front())}),
                "Relationship", "reaction-equilibrium", "final reaction gap is below one millith",
                "reference_or_one", 1.0e-3, 0.0, gap <= 1.0e-3,
                "s11-final-gap", "final", "\\Delta", "T,P,\\mathbf{x},A", "\\Delta\\leq10^{-3}", "pressure");
        }
    }
}

void check_s13() {
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    const auto problem = equilibrium_support::flash(2);
    const auto result = epcsaft::equilibrium::solve_equilibrium(model, problem);
    const double limit = result.convergence.criterion == nullptr
        ? 0.0 : convergence_limit(result.convergence);
    recording_support::compare(
        "phase-count", "phase_count", "count", static_cast<double>(result.phase_count), 1.0,
        recording_support::fields({recording_support::field("T", 300.0),
                                   recording_support::field("P", 1.0e5),
                                   recording_support::field("feed", problem.feed.values),
                                   recording_support::field("limit", limit)}),
        "Relationship", "single-phase", "single-phase solve returns one phase",
        "reference_or_one", 0.0, 0.0, result.phase_count == 1,
        "s13-phase-count", "phase", "N_{phase}", "T,P,feed", "N_{phase}=1", "phase");
    CHECK(result.success && result.phase_count == 1);
    CHECK(result.mole_fractions.size() == 2);
    recording_support::compare(
        "composition-0", "mole_fraction", "1", result.mole_fractions[0], problem.feed.values[0],
        recording_support::fields({recording_support::field("T", 300.0), recording_support::field("P", 1.0e5),
                                   recording_support::field("feed", problem.feed.values), recording_support::field("limit", limit),
                                   recording_support::field("component", 0U)}),
        "Relationship", "feed-balance", "phase composition follows feed", "reference_or_one", limit, 0.0,
        std::abs(result.mole_fractions[0] - problem.feed.values[0]) <= limit,
        "s13-composition", "component-0", "x_0", "T,P,feed,component", "x_0=x_{feed,0}", "component");
    recording_support::compare(
        "composition-1", "mole_fraction", "1", result.mole_fractions[1], problem.feed.values[1],
        recording_support::fields({recording_support::field("T", 300.0), recording_support::field("P", 1.0e5),
                                   recording_support::field("feed", problem.feed.values), recording_support::field("limit", limit),
                                   recording_support::field("component", 1U)}),
        "Relationship", "feed-balance", "phase composition follows feed", "reference_or_one", limit, 0.0,
        std::abs(result.mole_fractions[1] - problem.feed.values[1]) <= limit,
        "s13-composition", "component-1", "x_1", "T,P,feed,component", "x_1=x_{feed,1}", "component");
    CHECK(std::abs(result.mole_fractions[0] - problem.feed.values[0]) <= limit);
    CHECK(std::abs(result.mole_fractions[1] - problem.feed.values[1]) <= limit);
    CHECK(result.molar_densities.size() == 1 && result.molar_densities.front() > 0.0);
    const auto root = epcsaft::solve_density(
        model, 300.0, 1.0e5, problem.feed.values, epcsaft::DensityBranch::Liquid, std::nullopt);
    const double density_limit = std::max(limit, convergence_limit(root.convergence));
    const double specific_volume_difference =
        std::abs(1.0 / result.molar_densities.front() - 1.0 / root.molar_density);
    recording_support::compare(
        "density-root", "specific_volume_difference", "m3/mol", specific_volume_difference, 0.0,
        recording_support::fields({recording_support::field("T", 300.0),
                                   recording_support::field("P", 1.0e5),
                                   recording_support::field("feed", problem.feed.values),
                                   recording_support::field("density_limit", density_limit),
                                   recording_support::field("rho_result", result.molar_densities.front()),
                                   recording_support::field("rho_root", root.molar_density)}),
        "Relationship", "density-root", "equilibrium density matches direct density root",
        "reference_or_one", density_limit, 0.0, specific_volume_difference <= density_limit,
        "s13-density-root", "liquid", "\\Delta v", "T,P,feed,branch", "\\Delta v\\leq\\ell", "density");
}

void check_c3_availability() {
    // NaCl at 1 mol/kg with a water-only vapor: the vapor vanishes, so the solve fails with its
    // diagnostics retained and its state cleared; a neutral single-phase solve keeps its state.
    epcsaft::HelmholtzModel failure_model(properties_support::aqueous_electrolyte());
    auto failure_problem = equilibrium_support::flash(3);
    failure_problem.phases = {
        equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid),
        equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor),
    };
    failure_problem.phases[1].support = {0};
    failure_problem.temperature.value = 298.15;
    failure_problem.pressure.value = 1.0e5;
    failure_problem.feed.kind = epcsaft::equilibrium::FeedKind::Molalities;
    failure_problem.feed.values = {0.0, 1.0, 1.0};
    failure_problem.feed.solvent = {0};
    failure_problem.feed.solvent_mass_kg = 1.0;
    const auto compiled = epcsaft::equilibrium::compile_problem(failure_model, failure_problem);
    epcsaft::core::Capture deficit_capture;
    const auto failure = epcsaft::equilibrium::solve_equilibrium(failure_model, failure_problem);
    const bool retained = failure.convergence.criterion != nullptr
        && std::isfinite(failure.convergence.residual)
        && failure.rows.size() == failure.residuals.size()
        && std::all_of(failure.residuals.begin(), failure.residuals.end(),
                       [](double value) { return std::isfinite(value); })
        && !failure.message.empty();
    const bool cleared = !failure.success && !failure.has_state()
        && !failure.temperature && !failure.pressure
        && failure.amounts.empty() && failure.mole_fractions.empty()
        && failure.molar_densities.empty() && failure.coordinates.empty();
    const bool failure_reason = failure.message.find("phase V vanished") != std::string::npos;
    // The vanished phase's KKT-scaled stationarity rows vanish with its amount, so the
    // classified residual meets the association floor while the topology predicate fails.
    const bool stalled = failure.convergence.stop == epcsaft::core::Stop::stalled
        || failure.convergence.stop == epcsaft::core::Stop::converged_to_floor;
    // Michelsen deficit of the vanished phase from its chemical rows.
    double expected_deficit = 1.0;
    bool driving_row = false;
    const double limit = failure.convergence.criterion == nullptr
        ? 0.0 : convergence_limit(failure.convergence);
    for (std::size_t row = 0; row < failure.rows.size(); ++row)
        if (failure.rows[row].kind == epcsaft::equilibrium::RowKind::Chemical
            && failure.rows[row].phase == 1) {
            expected_deficit -= std::exp(-failure.residuals[row]);
            driving_row = driving_row || std::abs(failure.residuals[row]) > limit;
        }
    CHECK(driving_row);
    const auto& trace = deficit_capture.trace();
    const auto deficit = std::find_if(trace.events.rbegin(), trace.events.rend(), [](const auto& event) {
            return event.span == 0 && std::string(event.name) == "michelsen_deficit";
        });
    const auto distance = std::find_if(trace.events.rbegin(), trace.events.rend(), [](const auto& event) {
            return event.span == 0 && std::string(event.name) == "phase_distance";
        });
    CHECK(deficit != trace.events.rend() && deficit->values.size() == 1);
    CHECK(distance != trace.events.rend() && distance->values.size() == 1);
    recording_support::compare(
        "deficit", "michelsen_deficit", "1", deficit->values.front(), expected_deficit,
        recording_support::fields({recording_support::field("T", 298.15),
                                   recording_support::field("P", 1.0e5),
                                   recording_support::field("feed", failure_problem.feed.values),
                                   recording_support::field("expected_deficit", expected_deficit)}),
        "Relationship", "Michelsen-deficit", "trace deficit equals recomputed phase deficit",
        "reference_or_one", 0.0, 1e-12,
        std::abs(deficit->values.front() - expected_deficit) <= 1.0e-12 * std::max(1.0, std::abs(expected_deficit)),
        "c3-deficit", "trace", "D", "T,P,feed", "D=1-\\sum_i\\exp(-r_i)", "deficit");
    const double distance_reference = compiled.feed_composition[1] + compiled.feed_composition[2];
    recording_support::compare_custom(
        "distance", "phase_distance", "1", distance->values.front(), distance_reference,
        recording_support::fields({recording_support::field("T", 298.15),
                                   recording_support::field("P", 1.0e5),
                                   recording_support::field("feed", failure_problem.feed.values),
                                   recording_support::field("reference_distance", distance_reference)}),
        "Relationship", "boundary-failure", "phase distance reaches the vanished-phase feed amount",
        "isfinite(model) && isfinite(reference) && model >= reference",
        std::isfinite(distance->values.front()) && std::isfinite(distance_reference)
            && distance->values.front() >= distance_reference,
        "c3-distance", "trace", "d", "T,P,feed", "d\\geq n_{vanished}", "trace");

    epcsaft::HelmholtzModel success_model(properties_support::neutral_binary());
    auto success_problem = equilibrium_support::flash(2);
    success_problem.reactions.clear();
    success_problem.phases.resize(1);
    const auto success = epcsaft::equilibrium::solve_equilibrium(success_model, success_problem);
    const bool success_state = success.success && success.has_state()
        && success.convergence.criterion != nullptr
        && std::isfinite(success.convergence.residual);
    recording_support::compare(
        "c3-availability", "failed-result", "1",
        (retained && cleared && failure_reason && stalled) ? 1.0 : 0.0, 1.0,
        recording_support::fields({
            recording_support::field("solver_status", failure.solver_status),
            recording_support::field("stop", epcsaft::core::stop_name(failure.convergence.stop)),
            recording_support::field("residual", failure.convergence.residual),
            recording_support::field("row_count", failure.rows.size()),
            recording_support::field("residual_count", failure.residuals.size()),
            recording_support::field("message", failure.message),
            recording_support::field("T", failure_problem.temperature.value),
            recording_support::field("P", failure_problem.pressure.value),
            recording_support::field("feed", failure_problem.feed.values),
        }),
        "Relationship", "unavailable-result", "failed solve retains diagnostics and clears state",
        "absolute", 0.0, 0.0,
        retained && cleared && failure_reason && stalled, "native.equilibrium.c3", "C3",
        "failed-result", "availability", "T,P,feed", "has_state=0 with diagnostics retained");
    recording_support::compare(
        "c3-availability", "success-state", "1", success_state ? 1.0 : 0.0, 1.0,
        recording_support::fields({
            recording_support::field("solver_status", success.solver_status),
            recording_support::field("stop", epcsaft::core::stop_name(success.convergence.stop)),
            recording_support::field("residual", success.convergence.residual),
            recording_support::field("T", success_problem.temperature.value),
            recording_support::field("P", success_problem.pressure.value),
            recording_support::field("feed", success_problem.feed.values),
        }),
        "Relationship", "available-result", "successful solve exposes state diagnostics",
        "absolute", 0.0, 0.0, success_state, "native.equilibrium.c3", "C3",
        "success-state", "availability", "T,P,feed", "has_state=1 with convergence retained");
    auto symmetric = equilibrium_support::synthetic_reactive_problem();
    symmetric.phases[0].composition_guess.clear();
    symmetric.phases[1].composition_guess.clear();
    bool refused = false;
    try { (void)epcsaft::equilibrium::solve_equilibrium(
        epcsaft::HelmholtzModel(equilibrium_support::synthetic_reactive_ternary()), symmetric); }
    catch (const epcsaft::core::DeclarationError& error) {
        refused = std::string(error.what()).find("composition_guess") != std::string::npos;
    }
    CHECK(refused);
}

void check_s16() {
    epcsaft::HelmholtzModel model(equilibrium_support::synthetic_reactive_ternary());
    auto problem = equilibrium_support::synthetic_reactive_problem();
    const auto first = epcsaft::equilibrium::solve_equilibrium(model, problem);
    CHECK(first.success);
    for (const double fraction : {0.25, 0.5, 0.75}) {
        for (std::size_t i = 0; i < 3; ++i)
            problem.feed.values[i] = (1.0 - fraction) * first.mole_fractions[i]
                + fraction * first.mole_fractions[3 + i];
        const auto result = epcsaft::equilibrium::solve_equilibrium(model, problem, first);
        CHECK(result.success);
        if (!result.success) continue;
        const double limit = convergence_limit(result.convergence);
        for (std::size_t i = 0; i < 6; ++i) {
            const double difference = std::abs(
                result.mole_fractions[i] - first.mole_fractions[i]);
            const auto point = recording_support::number(fraction) + "/" + std::to_string(i);
            recording_support::compare(
                (std::string("continuation/") + point).c_str(),
                "mole_fraction_difference", "1", difference, 0.0,
                recording_support::fields({recording_support::field("T", problem.temperature.value),
                                           recording_support::field("P", problem.pressure.value),
                                           recording_support::field("feed", problem.feed.values),
                                           recording_support::field("fraction", fraction),
                                           recording_support::field("component", i),
                                           recording_support::field("limit", limit),
                                           recording_support::field("solution", result.mole_fractions)}),
                "Relationship", "continuation", "continuation composition stays within convergence limit",
                "reference_or_one", limit, 0.0, difference <= limit,
                "s16-composition", point.c_str(), "\\Delta x_i",
                "T,P,feed,\\alpha,component", "\\Delta x_i\\leq\\ell", "component");
        }
        for (std::size_t phase = 0; phase < 2; ++phase) {
            double affinity = 0.0;
            for (std::size_t row = 0; row < result.rows.size(); ++row)
                if (result.rows[row].kind == epcsaft::equilibrium::RowKind::Chemical
                    && result.rows[row].phase == phase)
                    affinity += problem.reactions[0].stoichiometry[result.rows[row].index]
                        * result.residuals[row];
            recording_support::compare(
                (std::string("affinity/") + recording_support::number(fraction) + "/"
                 + std::to_string(phase)).c_str(),
                "reaction_affinity", "1", affinity, 0.0,
                recording_support::fields({recording_support::field("T", problem.temperature.value),
                                           recording_support::field("P", problem.pressure.value),
                                           recording_support::field("feed", problem.feed.values),
                                           recording_support::field("fraction", fraction),
                                           recording_support::field("phase", phase),
                                           recording_support::field("limit", limit)}),
                "Relationship", "reaction-equilibrium", "continuation reaction affinity is zero",
                "reference_or_one", limit, 0.0, std::abs(affinity) <= limit,
                "s16-affinity", (std::to_string(fraction) + "-" + std::to_string(phase)).c_str(),
                "A", "T,P,feed,\\alpha,phase", "A=0", "fraction,phase");
        }
    }
}

// ------------------------------------------------------------------ E2 solves

namespace e2 {

using epcsaft::equilibrium::AmountSpec;
using epcsaft::equilibrium::FreeValue;
using epcsaft::equilibrium::PhaseKind;

bool record_bool(const char* id, const char* group, const char* claim, bool verdict,
                 const std::string& inputs) {
    recording_support::compare(id, "verdict", "bool", verdict ? 1.0 : 0.0, 1.0, inputs,
        "Relationship", "equilibrium", claim, "absolute", 0.0, 0.0, verdict, group, id,
        "verdict", "inputs", claim, "case");
    return verdict;
}

// The classifying residual: KKT-scaled stationarity rows for free-amount cells, unscaled
// rows for pinned cells (owner decision 2026-09-16); the unscaled rows stay diagnostics.
double largest_row(const epcsaft::equilibrium::EquilibriumResult& result) {
    return result.convergence.residual;
}

std::string solve_inputs(const epcsaft::equilibrium::EquilibriumProblem& problem,
                         const epcsaft::equilibrium::EquilibriumResult& result) {
    return recording_support::fields({
        recording_support::field("feed", problem.feed.values),
        recording_support::field("T", problem.temperature.value),
        recording_support::field("P", problem.pressure.value),
        recording_support::field("solver_status", result.solver_status),
        recording_support::field("stop", epcsaft::core::stop_name(result.convergence.stop)),
        recording_support::field("residual", result.convergence.residual),
        recording_support::field("x", result.mole_fractions),
        recording_support::field("amounts", result.amounts),
        recording_support::field("rho", result.molar_densities)});
}

bool not_run_record(const epcsaft::equilibrium::EquilibriumResult& result,
                    const epcsaft::equilibrium::CompiledSystem& compiled) {
    return !result.success && result.rows.size() == compiled.rows.size()
        && result.unknown_count == compiled.unknown_count()
        && result.jacobian_rank == compiled.jacobian_rank
        && result.generic_nullity == compiled.generic_nullity && result.generic_nullity > 0
        && result.residuals.empty() && result.convergence.stop == epcsaft::core::Stop::infeasible
        && std::isnan(result.convergence.residual) && result.convergence.criterion == nullptr
        && !result.requested_tolerance_met.has_value() && !result.convergence.floor.has_value()
        && !result.convergence.amplification.has_value() && !result.has_state()
        && result.coordinates.empty() && result.multipliers.empty();
}

}  // namespace e2

void check_e2_thermal() {
    using namespace e2;
    const double R = epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    epcsaft::HelmholtzModel pair(identical_ethane_pair());
    bool isomerisation_ok = true;
    // Qualification rows use the moderate offset; the dilute legacy offset is retained below as
    // a documented limitation of the unpinned Gibbs route (owner decision, 2026-09-16).
    struct Offset final { double a; bool qualified; };
    for (const Offset offset : {Offset{moderate_isomerisation_a(), true}, Offset{DILUTE_ISOMERISATION_A, false}}) {
    for (const double sign : offset.qualified ? std::vector<double>{1.0, -1.0} : std::vector<double>{1.0}) {
        for (const double temperature : {330.0, 350.0, 370.0}) {
            auto problem = equilibrium_support::flash(2);
            problem.phases = {equilibrium_support::phase("V", PhaseKind::Vapor)};
            problem.feed.values = {0.5, 0.5};
            problem.temperature = FreeValue::fixed(temperature);
            problem.pressure = FreeValue::fixed(2.0e5);
            problem.reactions = {{"A=B", {-sign, sign}, std::nullopt, std::nullopt, legacy_isomerisation(sign, offset.a)}};
            const auto result = epcsaft::equilibrium::solve_equilibrium(pair, problem);
            const double ln_k = offset.a + 1200.0 / temperature + 0.4 * std::log(temperature) + 0.002 * temperature;
            const double expected = std::exp(ln_k) / (1.0 + std::exp(ln_k));
            const double observed = result.success ? result.mole_fractions[1] : std::numeric_limits<double>::quiet_NaN();
            const auto ethane_state = pair.state_at_pressure(temperature, 2.0e5, {1.0, 0.0}, epcsaft::DensityBranch::Vapor);
            const double density_gap = result.success
                ? std::abs(result.molar_densities[0] / ethane_state.molar_density - 1.0)
                : std::numeric_limits<double>::quiet_NaN();
            const auto id = std::string("e2-2/isomerisation/") + (offset.qualified ? "" : "dilute-limit/")
                + (sign > 0.0 ? "forward/" : "reversed/") + recording_support::number(temperature);
            // Unpinned chemical rows are recovered from the KKT gradient N x_i mu_i, so the
            // dilute participant's row is bounded only by row_tolerance / (N x_A).
            const double dilute_bound = result.convergence.criterion->value / (1.0 * (1.0 - expected));
            const auto inputs = recording_support::fields({solve_inputs(problem, result),
                recording_support::field("a", offset.a), recording_support::field("ln_K", ln_k),
                recording_support::field("x_B_expected", expected),
                recording_support::field("rho_ethane", ethane_state.molar_density),
                recording_support::field("density_gap", density_gap),
                recording_support::field("dilute_row_bound", dilute_bound)});
            const bool verdict = result.success && std::abs(observed - expected) <= 1.0e-10
                && largest_row(result) <= result.convergence.criterion->value && density_gap <= 1.0e-9;
            recording_support::compare(id.c_str(), "mole_fraction", "1", observed, expected, inputs,
                "Relationship", "reaction-equilibrium", offset.qualified
                    ? "identical-species isomerisation follows x_B=K/(1+K) with the EOS density"
                    : "documented limitation: dilute unpinned participant row exceeds the criterion",
                "strict_absolute", 1.0e-10, 0.0, verdict, "e2-2-thermal", id.c_str(), "x_B", "T,P,K(T)", "x_B=K/(1+K)", "case");
            if (offset.qualified) {
                isomerisation_ok = isomerisation_ok && verdict;
            } else {
                const bool bounded = std::isfinite(result.convergence.residual)
                    && result.convergence.residual <= dilute_bound;
                recording_support::compare_custom((id + "/row-bound").c_str(), "row_residual", "1",
                    result.convergence.residual, dilute_bound, inputs, "Relationship", "conditioning",
                    "dilute recovered chemical row stays within row_tolerance/(N x_A)",
                    "isfinite(model) && model <= reference", bounded, "e2-2-thermal", id.c_str(),
                    "|r|", "T,P,K(T)", "|r|\\leq tol/(N x_A)", "case");
                // The documented limit is asserted; the unavailable x_B row above is not.
                isomerisation_ok = isomerisation_ok && bounded;
            }
        }
    }
    }
    // A -> 2B with matching species thermochemistry and an explicitly preconverted correlation.
    const double tref = 300.0;
    const epcsaft::ReactionLogPolynomial pressure_law{-10.0 / R, 3000.0 / R, 10.0 / R, 0.0, tref, 250.0, 450.0, "pressure-basis-P0"};
    const epcsaft::ReactionLogPolynomial density_law{
        pressure_law.a - std::log(R * tref * 1.0 / 1.0e5), pressure_law.b, pressure_law.c - 1.0, 0.0,
        tref, 250.0, 450.0, epcsaft::equilibrium::EOS_STANDARD_STATE_ID};
    const auto p3_pressure = epcsaft::evaluate_reaction_properties(pressure_law, 350.0);
    const auto p3_density = epcsaft::evaluate_reaction_properties(density_law, 350.0);
    const double cp_density = R * (2.0 * 350.0 * p3_density.derivative + 350.0 * 350.0 * p3_density.second_derivative);
    const auto p3_inputs = recording_support::fields({recording_support::field("T", 350.0),
        recording_support::field("a_rho", density_law.a), recording_support::field("b_rho", density_law.b),
        recording_support::field("c_rho", density_law.c)});
    recording_support::compare("e2-2/a-to-2b/ln_K_rho", "ln_K", "1", p3_density.ln_k, 3.5505788308151547, p3_inputs,
        "Relationship", "standard-state-conversion", "preconverted density-basis correlation",
        "strict_absolute", 1.0e-12, 0.0, std::abs(p3_density.ln_k - 3.5505788308151547) <= 1.0e-12,
        "e2-2-thermal", "ln_K_rho", "\\ln K_\\rho", "T", "\\ln K_\\rho=3.5505788308151547", "case");
    recording_support::compare("e2-2/a-to-2b/dh_rho", "reaction_enthalpy", "J/mol", p3_density.enthalpy, -2410.061916353634, p3_inputs,
        "Relationship", "standard-state-conversion", "density-basis reaction enthalpy", "strict_absolute", 1.0e-9, 0.0,
        std::abs(p3_density.enthalpy + 2410.061916353634) <= 1.0e-9, "e2-2-thermal", "dh_rho", "\\Delta h_\\rho", "T", "\\Delta h_\\rho", "case");
    recording_support::compare("e2-2/a-to-2b/dcp_rho", "reaction_heat_capacity", "J/(mol K)", cp_density, 1.6855373818467605, p3_inputs,
        "Relationship", "standard-state-conversion", "density-basis reaction heat capacity from P3 derivatives", "strict_absolute", 1.0e-12, 0.0,
        std::abs(cp_density - 1.6855373818467605) <= 1.0e-12, "e2-2-thermal", "dcp_rho", "\\Delta c_{p,\\rho}", "T", "\\Delta c_{p,\\rho}", "case");
    CHECK(std::abs(p3_density.ln_k - 3.5505788308151547) <= 1.0e-12);
    CHECK(std::abs(p3_density.enthalpy + 2410.061916353634) <= 1.0e-9 && std::abs(cp_density - 1.6855373818467605) <= 1.0e-12);
    epcsaft::HelmholtzModel thermal(identical_ethane_pair(), {}, species_record(30.0, 20.0));
    epcsaft::HelmholtzModel ethane(properties_support::ethane());
    const auto pure = ethane.state_at_pressure(350.0, 2.0e5, {1.0}, epcsaft::DensityBranch::Vapor);
    const double compressibility = 2.0e5 / (pure.molar_density * R * 350.0);
    const double phi = std::exp(pure.residual_chemical_potential_over_rt[0].value()) / compressibility;
    const double a = std::exp(p3_pressure.ln_k) * 1.0e5 / (phi * 2.0e5);
    const double y = 2.0 * a / (a + std::sqrt(a * a + 4.0 * a));
    const double total = 1.5 / (2.0 - y);
    auto doubling = equilibrium_support::flash(2);
    doubling.phases = {equilibrium_support::phase("V", PhaseKind::Vapor)};
    doubling.phases[0].composition_guess = {1.0 - y, y};
    doubling.feed.values = {0.5, 0.5};
    doubling.temperature = FreeValue::fixed(350.0);
    doubling.pressure = FreeValue::fixed(2.0e5);
    doubling.reactions = {{"A=2B", {-1.0, 2.0}, std::nullopt, std::nullopt, density_law}};
    const auto doubling_result = epcsaft::equilibrium::solve_equilibrium(thermal, doubling);
    const double composition_gap = doubling_result.success ? std::abs(doubling_result.mole_fractions[1] - y) : std::numeric_limits<double>::quiet_NaN();
    const double amount_gap = doubling_result.success
        ? std::max(std::abs(doubling_result.amounts[0] * doubling_result.mole_fractions[0] - total * (1.0 - y)),
                   std::abs(doubling_result.amounts[0] * doubling_result.mole_fractions[1] - total * y))
        : std::numeric_limits<double>::quiet_NaN();
    const bool doubling_ok = doubling_result.success && composition_gap <= 1.0e-8 && amount_gap <= 1.0e-8
        && largest_row(doubling_result) <= epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value;
    recording_support::compare("e2-2/a-to-2b/x_B", "mole_fraction", "1",
        doubling_result.success ? doubling_result.mole_fractions[1] : std::numeric_limits<double>::quiet_NaN(), y,
        recording_support::fields({solve_inputs(doubling, doubling_result), recording_support::field("Z", compressibility),
            recording_support::field("phi", phi), recording_support::field("A", a), recording_support::field("N_expected", total),
            recording_support::field("amount_gap", amount_gap)}),
        "Relationship", "reaction-equilibrium", "A=2B composition follows the independent pure-state relation y^2/(1-y)=K_P P0/(phi P)",
        "strict_absolute", 1.0e-8, 0.0, doubling_ok, "e2-2-thermal", "a-to-2b", "x_B", "T,P,K_P,\\phi", "y^2/(1-y)=A", "case");
    auto unconverted = doubling;
    unconverted.reactions[0].correlation = pressure_law;
    std::string mismatch;
    try { (void)epcsaft::equilibrium::compile_problem(thermal, unconverted); }
    catch (const std::invalid_argument& error) { mismatch = error.what(); }
    record_bool("e2-2/a-to-2b/reference-basis-mismatch", "e2-2-thermal", "reference-basis-mismatch",
        mismatch.find("reference-basis-mismatch") != std::string::npos,
        recording_support::fields({recording_support::field("identity", pressure_law.standard_state_id), recording_support::field("message", mismatch)}));
    epcsaft::HelmholtzModel mismatched(identical_ethane_pair(), {}, species_record(30.0, 21.0));
    std::string inconsistency;
    try { (void)epcsaft::equilibrium::compile_problem(mismatched, doubling); }
    catch (const std::invalid_argument& error) { inconsistency = error.what(); }
    record_bool("e2-2/a-to-2b/heat-capacity-inconsistency", "e2-2-thermal", "reference inconsistency",
        inconsistency.find("reference inconsistency") != std::string::npos,
        recording_support::fields({recording_support::field("cp_B", 21.0), recording_support::field("message", inconsistency)}));
    // Solved charged reaction with a manufactured consistent standard state at an admitted EOS state.
    const std::vector<double> x_star = {0.98, 0.01, 0.01};
    const std::vector<double> nu = {2.0, -1.0, -1.0};
    const std::vector<double> n_star = {1.96, 0.02, 0.02};
    for (const std::vector<std::size_t>& order : {std::vector<std::size_t>{0, 1, 2}, std::vector<std::size_t>{2, 0, 1}}) {
        epcsaft::HelmholtzModel salt(equilibrium_support::permute_parameters(properties_support::aqueous_electrolyte(), order));
        const auto x = equilibrium_support::permute_vector(x_star, order);
        const auto state = salt.state_at_pressure(298.15, 1.0e5, x, epcsaft::DensityBranch::Liquid);
        const double rt = R * 298.15;
        std::vector<double> gibbs(3);
        for (std::size_t i = 0; i < 3; ++i)
            gibbs[i] = -rt * (std::log(state.molar_density * x[i]) + state.residual_chemical_potential_over_rt[i].value());
        const auto stoichiometry = equilibrium_support::permute_vector(nu, order);
        double contraction = 0.0;
        for (std::size_t i = 0; i < 3; ++i) contraction += stoichiometry[i] * gibbs[i];
        auto problem = equilibrium_support::flash(3);
        problem.feed.values = equilibrium_support::permute_vector(std::vector<double>{0.0, 1.0, 1.0}, order);
        problem.temperature = FreeValue::fixed(298.15);
        problem.phases[0].support = {0, 1, 2};
        problem.phases[0].composition_guess = x;
        problem.reactions = {{"2W=Na+Cl", stoichiometry, std::exp(-contraction / rt), std::nullopt}};
        problem.standard_gibbs_j_per_mol = gibbs;
        const auto result = epcsaft::equilibrium::solve_equilibrium(salt, problem);
        const auto expected = equilibrium_support::permute_vector(n_star, order);
        double amount_gap = std::numeric_limits<double>::quiet_NaN(), charge = std::numeric_limits<double>::quiet_NaN();
        if (result.success) {
            amount_gap = 0.0; charge = 0.0;
            for (std::size_t i = 0; i < 3; ++i) {
                amount_gap = std::max(amount_gap, std::abs(result.amounts[0] * result.mole_fractions[i] - expected[i]));
                charge += salt.parameters().charge_number[i] * result.mole_fractions[i];
            }
        }
        // Rows use the suite's solved-row comparator max(criterion, floor): the association
        // gradient floor limits this ionic state (converged_to_floor), as in S1.
        const double limit = result.convergence.criterion == nullptr ? 0.0 : convergence_limit(result.convergence);
        const bool ok = result.success && amount_gap <= 1.0e-8 && std::abs(charge) <= 1.0e-10 && largest_row(result) <= limit;
        const auto id = std::string("e2-2/charged-solve/") + (order[0] == 0 ? "identity" : "permuted");
        recording_support::compare(id.c_str(), "amount_gap", "mol", ok ? amount_gap : std::numeric_limits<double>::quiet_NaN(), 0.0,
            recording_support::fields({solve_inputs(problem, result), recording_support::field("n_star", expected),
                recording_support::field("g_manufactured", gibbs), recording_support::field("K", std::exp(-contraction / rt)),
                recording_support::field("charge", charge), recording_support::field("row_limit", limit),
                recording_support::field("floor", epcsaft::core::optional_number(result.convergence.floor)),
                recording_support::field("order", std::vector<double>(order.begin(), order.end()))}),
            "Relationship", "reaction-equilibrium", "shared driver reproduces the manufactured charged solution n*",
            "strict_absolute", 1.0e-8, 0.0, ok, "e2-2-thermal", id.c_str(), "\\Delta n", "T,P,feed,g,K", "\\Delta n\\leq10^{-8}", "case");
    }
    CHECK(isomerisation_ok);
}

// Exact elimination and extensive homogeneity, checked against complete Helmholtz derivatives.
void check_pure_saturation() {
    using namespace epcsaft;
    using namespace epcsaft::equilibrium;
    const auto compare = [](const std::string& id, double actual, double reference, double scale) {
        const double error = std::abs(actual - reference) / std::max({scale, std::abs(actual), std::abs(reference)});
        recording_support::compare_custom(id.c_str(), "pure_saturation_identity", "scaled", actual, reference,
            recording_support::fields({recording_support::field("scale", scale), recording_support::field("error", error)}),
            "Relationship", "pure_isothermal_projection", "exact reduction, lift or homogeneity", "quantity scale",
            std::isfinite(error) && error <= 1e-10, "pure-saturation");
    };
    for (const auto& [name, parameters, temperature] : std::vector<std::tuple<std::string, Parameters, double>>{
        {"propane", caloric_specification_checks::propane_parameters(), 300},
        {"propane-low-pressure", caloric_specification_checks::propane_parameters(), 110},
        {"propane-packing-conversion", caloric_specification_checks::propane_parameters(), 216.4},
        {"acetone", pure_saturation_checks::acetone(), 300}, {"methanol", pure_saturation_checks::methanol(), 300}}) {
        HelmholtzModel model(parameters);
        const auto p = pure_saturation_checks::saturation(temperature);
        const auto c = compile_problem(model, p); const auto result = solve_compiled(c);
        const auto dedicated = solve_saturation(c);
        CHECK(dedicated.success && dedicated.requested_tolerance_met.value_or(false));
        CHECK(dedicated.residuals == result.residuals && dedicated.molar_densities == result.molar_densities);
        CHECK(dedicated.pressure == result.pressure && dedicated.phase_pressures.size() == 2);
        for (std::size_t i = 0; i < 2; ++i)
            CHECK(dedicated.phase_pressures[i] == model.pure_isothermal(temperature, 1 / dedicated.molar_densities[i]).pressure);
        if (!result.has_state()) throw std::runtime_error(name + ": " + result.message);
        for (double shift : {0.0, -0.01}) {
            const std::vector<double> u{std::log(result.molar_densities[0]) + shift, std::log(result.molar_densities[1]) - shift};
            const auto reduced = one_component_two_density_point(c, u); const auto& limited = reduced.phase;
            const auto repeated = one_component_two_density_point(c, u, nullptr, &limited,
                std::log(model.packing_fraction_coefficient(temperature, {1.0})));
            CHECK(repeated.phase.point == limited.point && repeated.rows == reduced.rows);
            const auto full = c.prepare(limited.point); const auto rows = c.evaluate_at(full).values;
            const auto jacobian = c.jacobian_at(full), limited_jacobian = c.jacobian_at(limited);
            const std::string id = name + "/" + std::to_string(shift);
            for (std::size_t phase = 0; phase < 2; ++phase) {
                for (double amount : {1.0, 100.0}) {
                    const double volume = amount * full.volume[phase];
                    const auto a = model.pure_isothermal(temperature, volume, amount).derivative;
                    const auto b = model.evaluate(temperature, volume, {amount}, DerivativeOrder::Hessian); const auto width = a.gradient.size();
                    const std::string key = id + "/" + std::to_string(phase) + "/" + std::to_string(amount);
                    compare(key + "/energy", a.value, b.value, amount);
                    for (std::size_t i : {1U, 2U}) {
                        const double si = i == 1 ? volume : amount;
                        compare(key + "/gradient" + std::to_string(i), a.gradient[i] * si, b.gradient[i] * si, amount);
                        for (std::size_t j : {1U, 2U}) {
                            const double sj = j == 1 ? volume : amount;
                            compare(key + "/hessian" + std::to_string(i) + std::to_string(j), a.hessian[i * width + j] * si * sj,
                                b.hessian[i * width + j] * si * sj, amount);
                        }
                    }
                    CHECK(a.gradient_floor >= 0 && std::isfinite(a.gradient_floor));
                    CHECK(a.site_fractions.size() == b.site_fractions.size());
                }
            }
            for (std::size_t i = 0; i < jacobian.size(); ++i) compare(id + "/full-J" + std::to_string(i), limited_jacobian[i], jacobian[i], 1);
            Eigen::MatrixXd lift = Eigen::MatrixXd::Zero(c.unknown_count(), 2);
            const auto offsets = c.phase_offsets();
            for (std::size_t phase = 0; phase < 2; ++phase) lift(offsets[phase] + 1 + (p.phases[phase].is_free() ? 1 : 0), phase) = 1;
            const auto& vapor = limited.pure[1]; const double kv = vapor.pressure_density_derivative / (GAS_CONSTANT_J_PER_MOL_K * temperature);
            lift(c.primal_size - 1, 1) = vapor.pressure_density_derivative / limited.volume[1]; lift(c.primal_size, 1) = kv;
            const auto n = static_cast<Eigen::Index>(c.unknown_count());
            const Eigen::MatrixXd projected = Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(jacobian.data(), n, n) * lift;
            for (std::size_t row = 0; row < c.rows.size(); ++row) {
                const auto kind = c.rows[row].kind; const auto phase = c.rows[row].phase;
                const int equation = phase == 0 && kind == RowKind::Pressure ? 0 : phase == 0 && kind == RowKind::Chemical ? 1 : -1;
                compare(id + "/row" + std::to_string(row), equation < 0 ? 0 : reduced.rows[equation], rows[row], 1);
                for (std::size_t col = 0; col < 2; ++col)
                    compare(id + "/J" + std::to_string(row) + std::to_string(col), equation < 0 ? 0 : reduced.jacobian[2 * equation + col], projected(row, col), 1);
            }
        }
        auto changed = p; std::reverse(changed.phases.begin(), changed.phases.end());
        changed.phases[0].amount = AmountSpec::free(); changed.phases[1].amount = AmountSpec::pinned(0);
        changed.phases[0].name = "gas"; changed.phases[1].name = "condensed"; changed.feed.values = {100};
        const auto lifted = solve_saturation(compile_problem(model, changed));
        if (!lifted.has_state()) throw std::runtime_error(name + " reversed: " + lifted.message);
        compare(name + "/feed-pressure", *lifted.pressure, *result.pressure, *result.pressure);
        compare(name + "/extensive-objective", lifted.objective, 100 * result.objective, 100);
        for (std::size_t phase = 0; phase < 2; ++phase) compare(name + "/phase-order" + std::to_string(phase), lifted.molar_densities[phase], result.molar_densities[1 - phase], result.molar_densities[1 - phase]);
        compare(name + "/free-amount", lifted.amounts[0], 100, 100);
    }
}

void check_e2_solve() {
    check_pure_saturation();
    using namespace e2;
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    auto one = equilibrium_support::flash(2);
    one.feed.values = {0.5, 0.5};
    one.phases[0].support = {0, 1};
    one.phases[0].amount = AmountSpec::pinned(1.0);
    one.temperature = FreeValue::variable(300.0);
    auto two = equilibrium_support::boundary();
    two.feed.values = {0.5, 0.5};
    two.phases[0].amount = AmountSpec::pinned(0.5);
    two.phases[1].amount = AmountSpec::pinned(0.5);
    two.temperature = FreeValue::variable(300.0);
    for (const auto* item : {&one, &two}) {
        const auto compiled = epcsaft::equilibrium::compile_problem(model, *item);
        const auto result = epcsaft::equilibrium::solve_equilibrium(model, *item);
        const auto id = std::string("e2-6/all-pinned-not-run/") + std::to_string(item->phases.size());
        record_bool(id.c_str(), "e2-6-solve", "nonreactive all-pinned declaration returns the exact pre-solve NOT_RUN record",
            compiled.generic_nullity > 0 && not_run_record(result, compiled),
            recording_support::fields({recording_support::field("phases", item->phases.size()),
                recording_support::field("rows", compiled.rows.size()), recording_support::field("unknowns", compiled.unknown_count()),
                recording_support::field("jacobian_rank", compiled.jacobian_rank), recording_support::field("generic_nullity", compiled.generic_nullity),
                recording_support::field("message", result.message)}));
    }
    // Positive pin in both orientations, then the count-deficient fixed-T/fixed-P form.
    epcsaft::HelmholtzModel ternary(properties_support::neutral_ternary());
    auto bubble = equilibrium_support::flash(3);
    bubble.phases = {equilibrium_support::phase("L", PhaseKind::Liquid),
                     equilibrium_support::phase("V", PhaseKind::Vapor, AmountSpec::pinned(0.0))};
    bubble.temperature.value = 180.0;
    bubble.pressure = FreeValue::variable(1.0e5);
    const auto bubble_result = epcsaft::equilibrium::solve_equilibrium(ternary, bubble);
    CHECK(bubble_result.success);
    for (const bool free_pressure : {true, false}) {
        auto pinned = bubble;
        pinned.phases[1].amount = AmountSpec::pinned(0.01);
        pinned.temperature = free_pressure ? FreeValue::fixed(180.0) : FreeValue::variable(180.0);
        pinned.pressure = free_pressure ? FreeValue::variable(*bubble_result.pressure)
                                        : FreeValue::fixed(*bubble_result.pressure);
        const auto compiled = epcsaft::equilibrium::compile_problem(ternary, pinned);
        const auto result = epcsaft::equilibrium::solve_equilibrium(ternary, pinned, bubble_result);
        const auto id = std::string("e2-6/positive-pin/") + (free_pressure ? "free-P" : "free-T");
        const bool ok = result.success && compiled.pinned_count == 1 && compiled.free_intensive_count == 1
            && largest_row(result) <= epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value
            && result.amounts[1] == 0.01 * compiled.feed_total;
        recording_support::compare(id.c_str(), "row_residual", "1", largest_row(result), 0.0,
            recording_support::fields({solve_inputs(pinned, result), recording_support::field("pin", 0.01)}),
            "Relationship", "equilibrium-residual", "one positive pin with one free intensive solves through the shared driver",
            "strict_absolute", epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value, 0.0, ok, "e2-6-solve", id.c_str(),
            "|r|", "T,P,feed,pin", "|r|\\leq10^{-10}", "case");
    }
    auto overdetermined = bubble;
    overdetermined.phases[1].amount = AmountSpec::pinned(0.01);
    overdetermined.pressure = FreeValue::fixed(*bubble_result.pressure);
    std::string message;
    try { (void)epcsaft::equilibrium::compile_problem(ternary, overdetermined); }
    catch (const std::invalid_argument& error) { message = error.what(); }
    record_bool("e2-6/positive-pin/fixed-T-fixed-P", "e2-6-solve", "positive pin with fixed T and P is a named structural rejection",
        message.find("pinned") != std::string::npos && message.find("free T or P") != std::string::npos,
        recording_support::fields({recording_support::field("message", message)}));
    // E2-7: explicit kind under a renamed phase, the single pure vapor phase, and the caller criterion.
    auto named = equilibrium_support::flash(2);
    auto renamed = named;
    renamed.phases[0].name = "liquid-phase";
    const auto named_result = epcsaft::equilibrium::solve_equilibrium(model, named);
    const auto renamed_result = epcsaft::equilibrium::solve_equilibrium(model, renamed);
    record_bool("e2-7/rename", "e2-7-criterion", "renaming an explicitly typed phase changes nothing physical",
        named_result.success && renamed_result.success && named_result.mole_fractions == renamed_result.mole_fractions
            && named_result.molar_densities == renamed_result.molar_densities && named_result.rows.size() == renamed_result.rows.size(),
        solve_inputs(renamed, renamed_result));
    epcsaft::HelmholtzModel ethane(properties_support::ethane());
    auto ethane_bubble = equilibrium_support::flash(1);
    ethane_bubble.phases = {equilibrium_support::phase("L", PhaseKind::Liquid),
                           equilibrium_support::phase("V", PhaseKind::Vapor, AmountSpec::pinned(0.0))};
    ethane_bubble.temperature.value = 108.6365;
    ethane_bubble.pressure = FreeValue::variable(1.0e5);
    const auto ethane_compiled = epcsaft::equilibrium::compile_problem(ethane, ethane_bubble);
    const auto ethane_result = epcsaft::equilibrium::solve_equilibrium(ethane, ethane_bubble);
    recording_support::compare("e2-6/ethane-bubble", "row_residual", "1", largest_row(ethane_result), 0.0,
        recording_support::fields({solve_inputs(ethane_bubble, ethane_result),
            recording_support::field("pressure_Pa", epcsaft::core::optional_number(ethane_result.pressure)),
            recording_support::field("jacobian_rank", ethane_compiled.jacobian_rank),
            recording_support::field("generic_nullity", ethane_compiled.generic_nullity)}),
        "Relationship", "equilibrium-residual", "the ethane bubble declaration has full rank and solves at 108.6365 K",
        "strict_absolute", epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value, 0.0,
        ethane_compiled.jacobian_rank == 7 && ethane_compiled.generic_nullity == 0 && ethane_result.success
            && ethane_result.has_state() && largest_row(ethane_result) <= epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value,
        "e2-6-solve", "ethane-bubble", "|r|", "T,feed,pins", "|r|\\leq10^{-10}", "case");
    auto vapor = equilibrium_support::flash(1);
    vapor.phases = {equilibrium_support::phase("V", PhaseKind::Vapor)};
    vapor.temperature.value = 320.0;
    const auto vapor_result = epcsaft::equilibrium::solve_equilibrium(ethane, vapor);
    const auto vapor_root = epcsaft::solve_density(ethane, 320.0, 1.0e5, {1.0}, epcsaft::DensityBranch::Vapor, std::nullopt);
    recording_support::compare("e2-7/pure-vapor", "row_residual", "1", largest_row(vapor_result), 0.0,
        recording_support::fields({solve_inputs(vapor, vapor_result), recording_support::field("root_density", vapor_root.molar_density),
            recording_support::field("stable_root_count", vapor_root.stable_root_count.value_or(0))}),
        "Relationship", "equilibrium-residual", "a single declared pure phase solves from its unique density root",
        "strict_absolute", epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value, 0.0,
        vapor_result.success && largest_row(vapor_result) <= epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value,
        "e2-7-criterion", "pure-vapor", "|r|", "T,P", "|r|\\leq10^{-10}", "case");
    std::size_t default_evaluations = 0;
    for (const double tolerance : {epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value, 1.0e-8, 1.0e-12, 1.0e-5}) {
        auto request = equilibrium_support::flash(2);
        request.feed.values = {0.5, 0.5};
        request.solver_options.values[epcsaft::policy::row] = tolerance;
        const auto result = epcsaft::equilibrium::solve_equilibrium(model, request);
        const double inner = tolerance * epcsaft::policy::EQUILIBRIUM_INNER_TOLERANCE_FACTOR.value;
        // Ipopt's numeric option storage round-trips decimal text (one ulp here).
        for (const auto* name : {"tol", "acceptable_tol", "constr_viol_tol", "dual_inf_tol"})
            CHECK(std::abs(result.solver_settings.at(name) - inner) <= 4 * std::numeric_limits<double>::epsilon() * inner);
        const bool default_request = tolerance == epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value;
        if (default_request) default_evaluations = result.evaluations;
        if (tolerance == 1.0e-5) CHECK(result.evaluations < default_evaluations && result.convergence.residual <= tolerance);
        const bool criterion_ok = result.convergence.criterion != nullptr
            && result.convergence.criterion->value == tolerance
            && (default_request ? result.convergence.criterion == &epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL
                                : std::string(result.convergence.criterion->name) == "EQUILIBRIUM_ROW_RESIDUAL");
        const bool flag_ok = std::isfinite(result.convergence.residual)
            ? (result.requested_tolerance_met.has_value()
               && *result.requested_tolerance_met == (result.convergence.residual <= tolerance))
            : !result.requested_tolerance_met.has_value();
        const bool warning_ok = !result.success || !result.requested_tolerance_met.has_value()
            || *result.requested_tolerance_met || result.message.find("requested row tolerance") != std::string::npos;
        const auto id = std::string("e2-7/tolerance/") + recording_support::number(tolerance);
        recording_support::compare(id.c_str(), "requested_tolerance_met", "bool",
            result.requested_tolerance_met.has_value() ? (*result.requested_tolerance_met ? 1.0 : 0.0) : std::numeric_limits<double>::quiet_NaN(),
            std::isfinite(result.convergence.residual) ? (result.convergence.residual <= tolerance ? 1.0 : 0.0) : std::numeric_limits<double>::quiet_NaN(),
            recording_support::fields({solve_inputs(request, result), recording_support::field("row_tolerance", tolerance), recording_support::field("inner_tol", result.solver_settings.at("tol")),
                recording_support::field("evaluations", result.evaluations), recording_support::field("iterations", result.iterations),
                recording_support::field("criterion", result.convergence.criterion == nullptr ? "absent" : result.convergence.criterion->name),
                recording_support::field("floor", epcsaft::core::optional_number(result.convergence.floor)),
                recording_support::field("success", result.success ? 1U : 0U), recording_support::field("message", result.message)}),
            "Relationship", "caller-criterion", "requested_tolerance_met is exactly recomputed_residual <= row_tolerance",
            "absolute", 0.0, 0.0, criterion_ok && flag_ok && warning_ok && (!default_request || result.success),
            "e2-7-criterion", id.c_str(), "met", "tolerance,residual", "met=(|r|\\leq tol)", "tolerance");
    }
    // Owned caller criterion survives copies through the trace after the source record is gone,
    // on the existing finite-floor fixture (aqueous electrolyte, association gradient floor).
    std::optional<epcsaft::core::Convergence> copied;
    {
        epcsaft::core::Capture capture;
        epcsaft::HelmholtzModel salt(properties_support::aqueous_electrolyte());
        auto request = equilibrium_support::flash(3);
        request.feed.values = {0.98, 0.01, 0.01};
        request.temperature.value = 298.15;
        request.solver_options.values[epcsaft::policy::row] = 1.0e-8;
        {
            const auto result = epcsaft::equilibrium::solve_equilibrium(salt, request);
            CHECK(result.success && result.convergence.floor.has_value());
        }
        for (const auto& span : capture.trace().spans)
            if (std::string(span.name) == "solve_equilibrium" && span.convergence.has_value()) copied = span.convergence;
    }
    const bool lifetime_ok = copied.has_value() && copied->criterion != nullptr && copied->criterion_owner
        && copied->floor.has_value() && std::isfinite(*copied->floor)
        && copied->criterion == copied->criterion_owner.get() && copied->criterion->value == 1.0e-8
        && std::string(copied->criterion->name) == "EQUILIBRIUM_ROW_RESIDUAL"
        && epcsaft::core::classify(copied->residual, *copied->criterion, copied->floor, copied->contraction, false,
                                   copied->amplification) == copied->stop;
    record_bool("e2-7/owned-criterion-lifetime", "e2-7-criterion",
        "a copied convergence record classifies identically after the producing result is destroyed", lifetime_ok,
        recording_support::fields({recording_support::field("residual", copied ? copied->residual : std::numeric_limits<double>::quiet_NaN()),
            recording_support::field("floor", copied ? epcsaft::core::optional_number(copied->floor) : "absent"),
            recording_support::field("stop", copied ? epcsaft::core::stop_name(copied->stop) : "absent")}));
}

void check_e2_reuse() {
    using namespace e2;
    struct Case final { const char* label; epcsaft::Parameters parameters; epcsaft::equilibrium::EquilibriumProblem problem; bool solved; };
    auto pinned_problem = equilibrium_support::boundary();
    pinned_problem.temperature.value = 180.0;
    auto free_t = pinned_problem;
    free_t.temperature = FreeValue::variable(180.0);
    free_t.pressure = FreeValue::fixed(1.0e5);
    auto reactive = free_t;
    reactive.reactions = {{"A-to-B", {-1.0, 1.0}, std::nullopt, 850.0}};
    auto associating = free_t;
    associating.temperature.value = 300.0;
    auto regular = equilibrium_support::flash(3);
    regular.phases = {equilibrium_support::phase("V", PhaseKind::Vapor), equilibrium_support::phase("L", PhaseKind::Liquid)};
    regular.temperature.value = 180.0;
    const std::vector<Case> cases = {
        {"regular", properties_support::neutral_ternary(), regular, true},
        {"pinned", properties_support::neutral_binary(), pinned_problem, false},
        {"pinned-free-T", properties_support::neutral_binary(), free_t, false},
        {"pinned-reactive-T", properties_support::neutral_binary(), reactive, false},
        {"pinned-association-T", properties_support::associating_binary(), associating, false},
    };
    const double criterion = 1.0e-9;
    for (const auto& item : cases) {
        const epcsaft::HelmholtzModel model(item.parameters);
        const auto compiled = epcsaft::equilibrium::compile_problem(model, item.problem);
        std::vector<double> a;
        if (item.solved) {
            const auto solved = epcsaft::equilibrium::solve_equilibrium(model, item.problem);
            CHECK(solved.success);
            a.assign(solved.coordinates.begin(), solved.coordinates.begin() + static_cast<std::ptrdiff_t>(compiled.primal_size));
        } else {
            a = compiled.generic_point();
            a[compiled.primal_size - compiled.free_intensive_count] = item.problem.temperature.free ? 193.0 : a[compiled.primal_size - compiled.free_intensive_count];
        }
        epcsaft::equilibrium::GibbsMinimizer nlp(compiled, a);
        Ipopt::Index variables = 0, constraints = 0, jacobian_entries = 0, hessian_entries = 0;
        Ipopt::TNLP::IndexStyleEnum index_style;
        CHECK(nlp.get_nlp_info(variables, constraints, jacobian_entries, hessian_entries, index_style));
        const auto n = variables;
        const auto m = constraints;
        CHECK(static_cast<std::size_t>(n) == a.size());
        std::vector<double> lambda_one(static_cast<std::size_t>(m)), lambda_two(static_cast<std::size_t>(m));
        for (std::size_t row = 0; row < lambda_one.size(); ++row) {
            lambda_one[row] = 0.35 + 0.1 * static_cast<double>(row);
            lambda_two[row] = -0.2 + 0.07 * static_cast<double>(row);
        }
        auto b = a;
        for (double& value : b) value += 1.0e-3;
        const auto unit = callback_unit(compiled, a);  // fixed by the start a for every point below
        std::vector<Ipopt::Index> rows_index(static_cast<std::size_t>(hessian_entries)), cols_index(static_cast<std::size_t>(hessian_entries));
        CHECK(nlp.eval_h(n, a.data(), true, 1.0, m, lambda_one.data(), true, hessian_entries, rows_index.data(), cols_index.data(), nullptr));
        const auto callback = [&](std::vector<double> point, double objective_factor, const std::vector<double>& lambda) {
            for (Ipopt::Index i = 0; i < n; ++i) point[i] /= unit[i];
            Ipopt::Number f = 0.0;
            std::vector<Ipopt::Number> g(static_cast<std::size_t>(m)), grad(static_cast<std::size_t>(n));
            std::vector<Ipopt::Number> jac(static_cast<std::size_t>(jacobian_entries)), hess(static_cast<std::size_t>(hessian_entries));
            CHECK(nlp.eval_f(n, point.data(), true, f));
            CHECK(nlp.eval_g(n, point.data(), false, m, g.data()));
            CHECK(nlp.eval_grad_f(n, point.data(), false, grad.data()));
            CHECK(nlp.eval_jac_g(n, point.data(), false, m, jacobian_entries, nullptr, nullptr, jac.data()));
            CHECK(nlp.eval_h(n, point.data(), false, objective_factor, m, lambda.data(), true, hessian_entries, nullptr, nullptr, hess.data()));
            for (Ipopt::Index i = 0, entry = 0; i < n; ++i) {  // restore S grad, J S and S H S to physical coordinates
                grad[i] /= unit[i];
                for (Ipopt::Index row = 0; row < m; ++row) jac[row * n + i] /= unit[i];
                for (Ipopt::Index column = 0; column <= i; ++column) hess[entry++] /= unit[i] * unit[column];
            }
            std::vector<double> flat = {f};
            flat.insert(flat.end(), g.begin(), g.end());
            flat.insert(flat.end(), grad.begin(), grad.end());
            flat.insert(flat.end(), jac.begin(), jac.end());
            flat.insert(flat.end(), hess.begin(), hess.end());
            return flat;
        };
        const auto fresh = [&](const std::vector<double>& point, double objective_factor, const std::vector<double>& lambda) {
            const auto evaluation = compiled.evaluate(point);
            const auto gradient = compiled.objective_gradient(point);
            const auto jacobian = compiled.jacobian(point);
            const auto objective = compiled.hessian(point, {});
            const auto weighted = compiled.hessian(point, lambda);
            const auto active = callback_rows(compiled);
            std::vector<double> flat = {compiled.pinned_count == 0 ? evaluation.objective : 0.0};
            for (const std::size_t row : active) flat.push_back(evaluation.values[row]);
            flat.insert(flat.end(), gradient.begin(), gradient.end());
            for (const std::size_t row : active)
                for (Ipopt::Index column = 0; column < n; ++column) flat.push_back(jacobian[row * static_cast<std::size_t>(n) + static_cast<std::size_t>(column)]);
            for (Ipopt::Index row = 0; row < n; ++row)
                for (Ipopt::Index column = 0; column <= row; ++column) {
                    const std::size_t index = static_cast<std::size_t>(row) * static_cast<std::size_t>(n) + static_cast<std::size_t>(column);
                    flat.push_back(objective_factor * objective[index] + (weighted[index] - objective[index]));
                }
            return flat;
        };
        const auto compare_flat = [&](const char* step, const std::vector<double>& left, const std::vector<double>& right) {
            CHECK(left.size() == right.size());
            double worst = 0.0;
            for (std::size_t i = 0; i < left.size(); ++i)
                worst = std::max(worst, std::abs(left[i] - right[i]) / std::max(1.0, std::abs(right[i])));
            const auto id = std::string("e2-reuse/") + item.label + "/" + step;
            recording_support::compare(id.c_str(), "callback_difference", "1", worst, 0.0,
                recording_support::fields({recording_support::field("case", item.label), recording_support::field("step", step),
                    recording_support::field("entries", left.size()), recording_support::field("point", a)}),
                "Relationship", "current-point-reuse", "cached callback values equal fresh evaluation at the same point",
                "reference_or_one", criterion, 0.0, worst <= criterion, "e2-reuse", id.c_str(), "\\Delta", "z,\\lambda,\\sigma",
                "\\Delta\\leq10^{-9}", "case,step");
        };
        const auto first = callback(a, 0.5, lambda_one);
        compare_flat("a-first", first, fresh(a, 0.5, lambda_one));
        compare_flat("a-changed-multipliers", callback(a, 1.0, lambda_two), fresh(a, 1.0, lambda_two));
        compare_flat("a-zero-objective", callback(a, 0.0, lambda_one), fresh(a, 0.0, lambda_one));
        compare_flat("b", callback(b, 0.5, lambda_one), fresh(b, 0.5, lambda_one));
        compare_flat("a-return", callback(a, 0.5, lambda_one), first);
    }
}

}  // namespace equilibrium_checks

namespace equilibrium_checks {
void check_ideal_gas() {
    using namespace epcsaft;
    using namespace epcsaft::equilibrium;
    Parameters parameters;
    parameters.ideal_gas = true;
    parameters.molar_mass_kg_per_mol = {0.02, 0.03};
    parameters.charge_number = {0, 0};
    ThermochemistryRecord record{100000.0, {}};
    for (double cp : {30.0, 40.0}) record.component_intervals.push_back({{400, 1000, true, true,
        {IdealPolynomial{{cp}, 600, 1000, 10}, 100000}}});
    HelmholtzModel model(parameters, {}, record);
    const double T = 600, P = 100000, rt = GAS_CONSTANT_J_PER_MOL_K * T;
    const auto state = model.state(T, P / rt, {0.25, 0.75});
    CHECK(close_with_tolerance(state.pressure / P, 1, 1e-10));
    CHECK(close_with_tolerance(state.pressure_density_derivative / rt, 1, 1e-10));
    const auto residual = model.evaluate(T, rt / P, {0.25, 0.75}, DerivativeOrder::Hessian);
    CHECK(residual.value == 0);
    for (auto value : residual.gradient) CHECK(value == 0);
    for (auto value : residual.hessian) CHECK(value == 0);
    auto problem = equilibrium_support::flash(2);
    problem.phases = {equilibrium_support::phase("V", PhaseKind::Vapor)};
    problem.temperature = FreeValue::fixed(T); problem.pressure = FreeValue::fixed(P);
    problem.feed.values = {0.25, 0.75};
    problem.solver_options.values[policy::row] = 1e-10;
    const auto c = compile_problem(model, problem);
    auto q = c.generic_point();
    const auto layout = c.prepare(q).shape.front();
    q[layout.x[0]] = std::log(0.2); q[layout.x[1]] = std::log(0.6);
    CHECK(layout.amount.has_value());
    q[*layout.amount] = 2;
    const double v = 2 * rt / P;
    q[layout.volume] = -std::log(v);
    // No reactions or supplied standards: the compiler uses zero species standards even with thermal records.
    const auto actual = c.evaluate(q);
    const auto gradient = c.objective_gradient(q), hessian = c.hessian(q);
    const auto n = q.size();
    const auto jacobian = c.jacobian(q);
    bool chemical_identity = true, chemical_difference = true;
    for (std::size_t j = 0; j < 2; ++j) {
        auto plus = q, minus = q;
        constexpr double step = 1e-5;
        plus[layout.x[j]] += step; minus[layout.x[j]] -= step;
        const auto up = c.evaluate(plus).values, down = c.evaluate(minus).values;
        for (std::size_t row = 0; row < c.rows.size(); ++row) if (c.rows[row].kind == RowKind::Chemical) {
            const double value = jacobian[row * n + layout.x[j]];
            const double identity = c.rows[row].index == j ? 1.0 : 0.0;
            const double difference = (up[row] - down[row]) / (2 * step);
            chemical_identity &= close_with_tolerance(value, identity, 1e-10);
            chemical_difference &= close_with_tolerance(value, difference, 1e-10);
            std::printf("coordinate,chemical-jacobian,i=%zu,j=%zu,actual=%.17g,reference=%.17g,central=%.17g\n",
                c.rows[row].index, j, value, identity, difference);
        }
    }
    const double expected = 2 * (0.2 * (std::log(0.2 / v) - 1) + 0.6 * (std::log(0.6 / v) - 1) + 2);
    CHECK(close_with_tolerance(actual.objective, expected, 1e-10));
    CHECK(close_with_tolerance(gradient[layout.volume], -2.4, 1e-10));
    CHECK(close_with_tolerance(hessian[layout.volume * n + layout.volume], 4, 1e-10));
    for (std::size_t i = 0; i < 2; ++i) {
        const double x = i == 0 ? 0.2 : 0.6;
        CHECK(close_with_tolerance(gradient[layout.x[i]], 2 * x * std::log(x / v), 1e-10));
        CHECK(close_with_tolerance(hessian[layout.x[i] * n + layout.volume], 2 * x, 1e-10));
        CHECK(close_with_tolerance(hessian[layout.x[i] * n + layout.x[i]], 2 * x * (std::log(x / v) + 1), 1e-10));
        std::printf("coordinate,composition-gradient,i=%zu,actual=%.17g,reference=%.17g,error=%.17g\n",i,gradient[layout.x[i]],2*x*std::log(x/v),gradient[layout.x[i]]-2*x*std::log(x/v));
        std::printf("coordinate,mixed-composition-volume,i=%zu,actual=%.17g,reference=%.17g,error=%.17g\n",i,hessian[layout.x[i]*n+layout.volume],2*x,hessian[layout.x[i]*n+layout.volume]-2*x);
        std::printf("coordinate,composition-curvature,i=%zu,actual=%.17g,reference=%.17g,error=%.17g\n",i,hessian[layout.x[i]*n+layout.x[i]],2*x*(std::log(x/v)+1),hessian[layout.x[i]*n+layout.x[i]]-2*x*(std::log(x/v)+1));
    }
    std::printf("coordinate,objective,actual=%.17g,reference=%.17g,error=%.17g\n",actual.objective,expected,actual.objective-expected);
    std::printf("coordinate,volume-gradient,actual=%.17g,reference=-2.4,error=%.17g\n",gradient[layout.volume],gradient[layout.volume]+2.4);
    std::printf("coordinate,volume-curvature,actual=%.17g,reference=4,error=%.17g\n",hessian[layout.volume*n+layout.volume],hessian[layout.volume*n+layout.volume]-4);
    const auto solved = solve_equilibrium(model, problem);
    std::printf("ideal-gas,success=%d,cause=%s,P=%.17g,slope=%.17g,F=%.17g,q_gradient=%.17g,qq=%.17g\n",
        solved.success, solved.message.c_str(), state.pressure, state.pressure_density_derivative,
        actual.objective, gradient[layout.volume], hessian[layout.volume*n+layout.volume]);
    CHECK(solved.success);
    CHECK(solved.amounts.size() == 1);
    CHECK(close_with_tolerance(solved.amounts[0], 1, 1e-10));
    for (std::size_t i = 0; i < 2; ++i) CHECK(std::abs(solved.amounts[0] * solved.mole_fractions[i] - problem.feed.values[i]) <= 2e-9);
    CHECK(!solved.residuals.empty());
    for (auto row : solved.residuals) CHECK(std::abs(row) <= 1e-10);
    auto target = problem;
    target.pressure = FreeValue::fixed(101000);
    const auto direct = solve_equilibrium(model, target);
    const auto path = continue_equilibrium(model, {problem, target, solved, {0.25, 1e-6, 20, 4}, {}});
    std::printf("ideal-continuation,success=%d,stop=%s,cause=%s,attempts=%zu,retries=%zu,accepted=%zu,reached=%.17g,message=%s\n",
        path.success, continuation_stop_name(path.stop), path.cause_stop ? continuation_stop_name(*path.cause_stop) : "none",
        path.attempts, path.retries, path.accepted_steps, path.reached_path, path.message.c_str());
    CHECK(chemical_identity); CHECK(chemical_difference);
    CHECK(direct.success && path.success && path.stop == ContinuationStop::reached_target);
    CHECK(path.reached_path == 1 && path.states.size() > 1 && path.states.size() == path.path_values.size());
    for (std::size_t k = 0; k < path.states.size(); ++k) {
        const auto& point = path.states[k];
        const double pressure = P + 1000 * path.path_values[k];
        std::printf("ideal-continuation,path=%.17g,P=%.17g,rho=%.17g,reference=%.17g\n",
            path.path_values[k], pressure, point.molar_densities[0], pressure / rt);
        CHECK(close_with_tolerance(point.molar_densities[0], pressure / rt, 1e-10));
        CHECK(!point.residuals.empty());
        for (double row : point.residuals) CHECK(std::abs(row) <= 1e-10);
        for (std::size_t i = 0; i < 2; ++i)
            CHECK(std::abs(point.amounts[0] * point.mole_fractions[i] - problem.feed.values[i]) <= 2e-9);
    }
    CHECK(close_with_tolerance(path.states.back().molar_densities[0], direct.molar_densities[0], 1e-10));
    bool refused = false;
    try { stability_test(c, solved, {}); } catch (const std::invalid_argument& e) { refused = std::string(e.what()).find("ideal-gas") != std::string::npos; std::printf("refusal,stability,T=600,P=100000,cause=%s\n",e.what()); }
    CHECK(refused);
    problem.phases.clear(); problem.topology_declared = false; problem.discovery = StabilityOptions{};
    for (bool direct : {false, true}) {
        refused = false;
        try { if (direct) discover_phases(model, problem); else solve_equilibrium(model, problem); }
        catch (const std::invalid_argument& e) { refused = std::string(e.what()).find("ideal-gas") != std::string::npos; std::printf("refusal,discovery,direct=%d,T=600,P=100000,cause=%s\n",direct,e.what()); }
        CHECK(refused);
    }
}
} // namespace equilibrium_checks

namespace equilibrium_checks {
void check_species_thermochemistry() {
    using namespace epcsaft;
    using namespace epcsaft::equilibrium;
    Parameters p; p.ideal_gas = true; p.molar_mass_kg_per_mol = {0.092011, 0.0460055}; p.charge_number = {0, 0};
    // NIST WebBook C10544726/C10102440, Gas Phase Heat Capacity (Shomate Equation), A-H rows.
    ThermochemistryRecord r{100000, {{{500,1000,true,true, {IdealShomate{{34.05274,191.9845,-151.0575,44.39350,-0.158949,-8.893428,293.7724,9.078988},9078.988},100000}}},
        {{298,1200,true,true, {IdealShomate{{16.10857,75.89525,-54.38740,14.30777,0.239423,26.17464,240.5386,33.09502},33095.02},100000}}}}};
    const std::vector<double> nu = {-1,2};
    auto problem = equilibrium_support::flash(2); problem.phases = {equilibrium_support::phase("V",PhaseKind::Vapor)};
    problem.temperature = FreeValue::fixed(500); problem.pressure = FreeValue::fixed(100000); problem.feed.values = {1,0};
    problem.reactions = {{"N2O4=2NO2",nu,{},{},{},{}}}; problem.solver_options.values[policy::row] = 1e-10;
    HelmholtzModel model(p,{},r); const auto c = compile_problem(model,problem); const auto solved = solve_equilibrium(model,problem);
    CHECK(solved.success); CHECK(!solved.residuals.empty());
    for(auto v:solved.residuals) CHECK(std::abs(v)<=1e-10);
    CHECK(std::abs(solved.amounts[0]*solved.mole_fractions[0]-0.0013392705057121423)<=2e-9);
    CHECK(std::abs(solved.amounts[0]*solved.mole_fractions[1]-1.9973214589885757)<=2e-9);
    const auto source = reaction_thermochemistry_density(r,nu,500);
    CHECK(source.status==ThermochemistryStatus::Available);
    CHECK(std::abs(source.pressure_gibbs+30375.9201906513)<=1e-6);
    CHECK(source.k_pressure && std::abs(*source.k_pressure/1490.3509615073785-1)<=1e-10);
    std::printf("species-target,T=500,G_P=%.17g,K_P=%.17g,G_rho=%.17g,n_N2O4=%.17g,n_NO2=%.17g,cause=%s\n",
        source.pressure_gibbs,*source.k_pressure,source.gibbs,solved.amounts[0]*solved.mole_fractions[0],solved.amounts[0]*solved.mole_fractions[1],solved.message.c_str());
    const auto read = reaction_activities(model, solved)[0];
    CHECK(read.log_quotient && std::abs(*read.log_quotient - source.ln_k_density) <= 1e-10);
    const double pressure_q = *read.log_quotient + std::log(GAS_CONSTANT_J_PER_MOL_K * 500 / 100000);
    CHECK(std::abs(pressure_q - source.ln_k_pressure) <= 1e-10);
    for (const auto& v : model.state(500, solved.molar_densities[0], solved.mole_fractions).residual_chemical_potential_over_rt) CHECK(v && *v == 0);
    auto candidate = solved;
    const auto snapshot_point = c.prepare(solved.coordinates);
    candidate.coordinates[snapshot_point.shape[0].volume] += 0.01;
    candidate.molar_densities[0] *= std::exp(0.01);
    const auto shifted = reaction_activities(model, candidate)[0];
    CHECK(shifted.log_quotient && std::abs(*shifted.log_quotient - *read.log_quotient - 0.01) <= 1e-10);
    for (std::size_t i = 0; i < 2; ++i)
        CHECK(std::abs(*shifted.species.log_activities[i] - *read.species.log_activities[i] - 0.01) <= 1e-10);
    check_activity_defect(model, c, candidate, source.ln_k_density);
    candidate = solved;
    const double extent = 1e-4, total = solved.amounts[0] + extent;
    for (std::size_t i = 0; i < 2; ++i) {
        candidate.mole_fractions[i] = (solved.amounts[0] * solved.mole_fractions[i] + nu[i] * extent) / total;
        candidate.coordinates[snapshot_point.shape[0].x[i]] = std::log(candidate.mole_fractions[i]);
    }
    candidate.amounts[0] = total; candidate.coordinates[*snapshot_point.shape[0].amount] = total;
    check_activity_defect(model, c, candidate, source.ln_k_density);
    std::printf("activity-gas,lnQ=%.17g,lnK=%.17g,pressureQ=%.17g,pressureK=%.17g,densityShift=%.17g\n",
        *read.log_quotient, source.ln_k_density, pressure_q, source.ln_k_pressure, *shifted.log_quotient - *read.log_quotient);
    for(bool missing:{false,true}) {
        auto unavailable = r; if(missing) unavailable.component_intervals[0].clear();
        const double T = missing?500:499;
        const auto value = reaction_thermochemistry_density(unavailable,nu,T);
        CHECK(value.status==(missing?ThermochemistryStatus::MissingRecord:ThermochemistryStatus::OutsideInterval));
        CHECK(!value.cause.empty()); CHECK(!value.k_pressure);
        HelmholtzModel absent(p,{},unavailable); auto request = problem; request.temperature = FreeValue::fixed(T);
        bool refused=false;try{compile_problem(absent,request);}catch(const std::invalid_argument&e){const std::string cause=e.what();refused=cause.find(value.cause)!=std::string::npos&&cause.find("N2O4=2NO2")!=std::string::npos;}
        CHECK(refused); std::printf("unavailable,%s\n",value.cause.c_str());
    }
    const auto fourth=reference_calorics(c,500,100000,true,false,false,4), caloric_source=reference_calorics(c,500,100000,true,false,false,2);
    CHECK(fourth.status==SolvedStateActionStatus::EnergyOrderUnavailable);
    CHECK(caloric_source.status==SolvedStateActionStatus::ReferenceUnavailable);
    std::printf("refusal,reference-calorics,order=4,T=500,P=100000,cause=%s\n",fourth.message.c_str());
    std::printf("refusal,reference-calorics,order=2,T=500,P=100000,cause=%s\n",caloric_source.message.c_str());
    for(auto kind:{SolvedStateObservableKind::SystemEnthalpy,SolvedStateObservableKind::EquilibriumIsobaricHeatCapacity,SolvedStateObservableKind::PhaseGibbsMolarVolume}) {
        SolvedStateActionRequest request{{{kind}},{{1,0,{0,0},{}}}};
        const auto action = solved_state_actions(c,solved,request);
        CHECK(action.results.size()==1&&action.results[0].status==SolvedStateActionStatus::EnergyOrderUnavailable);
        std::printf("refusal,caloric,kind=%d,T=500,P=100000,cause=%s\n",int(kind),action.results[0].diagnostic_message.value_or("").c_str());
    }
    for(std::size_t family:{std::size_t(0),std::size_t(2)}) for(std::size_t kind=0;kind<5;++kind) {
        SolvedStateObservable o; o.kind=SolvedStateObservableKind::PhaseProperty; o.property=property_observable(family,kind);
        const SolvedStateActionRequest request{{o},{{1,0,{0,0},{}}}};
        const auto action=solved_state_actions(c,solved,request);
        CHECK(action.results[0].status==SolvedStateActionStatus::EnergyOrderUnavailable);
        std::printf("refusal,phase-property,kind=%zu,T=500,P=100000,cause=%s\n",5*family+kind,action.results[0].diagnostic_message.value_or("").c_str());
    }
    for(bool trial:{false,true}) {
        std::vector<std::vector<double>> basis;
        bool refused=false;try{if(trial)compile_stability_trial(model,500,100000,problem.phases[0],c.n_feed,c.standard_gibbs,basis);else cooling_onset(c);}
        catch(const std::invalid_argument&e){refused=std::string(e.what()).find("ideal-gas")!=std::string::npos;std::printf("refusal,phase-search,trial=%d,T=500,P=100000,cause=%s\n",trial,e.what());} CHECK(refused);
    }
    for(std::size_t excluded=0;excluded<11;++excluded) {
        auto declaration=problem;
        if(excluded==0) declaration.temperature=FreeValue::variable(500);
        if(excluded==1) declaration.pressure=FreeValue::variable(100000);
        if(excluded==2) declaration.enthalpy_j=0;
        if(excluded==3) declaration.entropy_j_per_k=0;
        if(excluded==4) declaration.solids.push_back({});
        if(excluded==5) declaration.source_reference=NeutralReferenceDeclaration{};
        if(excluded==6) declaration.critical=CriticalDeclaration{};
        if(excluded==7) declaration.envelope=EnvelopeCondition{};
        if(excluded==8) declaration.azeotrope=true;
        if(excluded==9) declaration.phases.push_back(declaration.phases[0]);
        if(excluded==10) declaration.phases[0].kind=PhaseKind::Liquid;
        bool refused=false;try{compile_problem(model,declaration);}catch(const std::invalid_argument&e){refused=std::string(e.what()).find("ideal-gas")!=std::string::npos;std::printf("refusal,compiler,selection=%zu,T=500,P=100000,cause=%s\n",excluded,e.what());}
        CHECK(refused);
    }
    HelmholtzModel nonideal(properties_support::neutral_binary(),{},r);
    const auto nonideal_compiled=compile_problem(nonideal,problem);
    CHECK(reference_calorics(nonideal_compiled,500,100000,true,false,false,2).status==SolvedStateActionStatus::ReferenceUnavailable);
    const auto nonideal_solved=solve_equilibrium(nonideal,problem); CHECK(nonideal_solved.success);
    const auto deferred=solved_state_actions(nonideal_compiled,nonideal_solved,
        {{{SolvedStateObservableKind::EquilibriumIsobaricHeatCapacity}},{{1,0,{0,0},{}}}});
    CHECK(deferred.results[0].status==SolvedStateActionStatus::ReferenceUnavailable);
    CHECK(deferred.diagnostics.state_factorizations==0&&deferred.diagnostics.maximum_energy_order==0);
    std::printf("refusal,nonideal-species-caloric,T=500,P=100000,cause=%s\n",deferred.results[0].diagnostic_message.value_or("").c_str());
    for(bool entropy:{false,true}) {
        auto caloric=problem;if(entropy)caloric.entropy_j_per_k=0;else caloric.enthalpy_j=0;
        bool refused=false;try{compile_problem(nonideal,caloric);}catch(const std::invalid_argument&e){refused=std::string(e.what()).find("species-derived")!=std::string::npos;}CHECK(refused);
    }
    // Manufactured linear-Cp reference, independently integrated and differentiated.
    ThermochemistryRecord linear{100000,{}};
    linear.component_intervals={{{400,1000,true,true,{IdealPolynomial{{30,0.01},600,1000,10},100000}}},
                               {{400,1000,true,true,{IdealPolynomial{{40,0.03},600,5000,30},100000}}}};
    const double T=750,R=GAS_CONSTANT_J_PER_MOL_K,dt=T-600,L=std::log(R*T/100000);
    const double h=9000+50*dt+0.025*dt*dt,s=50+20*std::log(T/600)+0.05*dt,cp=50+0.05*dt;
    const std::array<double,4> expected={h-T*s+R*T*L,-s+R*(L+1),-cp/T+R/T,cp/(T*T)-0.05/T-R/(T*T)};
    const auto actual=reaction_thermochemistry_density(linear,nu,T);
    const auto targets=reaction_target_series(T,problem.reactions,3,&linear);
    HelmholtzModel manufactured(p,{},linear);auto request=problem;request.temperature=FreeValue::fixed(T);
    const auto compiled=compile_problem(manufactured,request); const auto point=compiled.generic_point(); const auto prepared=compiled.prepare(point);
    std::vector<std::vector<double>> fixed(4,std::vector<double>(point.size(),0.0)); fixed[0]=point;
    std::vector<double> thermal_control(compiled.control_width(),0.0); thermal_control[0]=1;
    const auto reduced=compiled.row_series(prepared,fixed,thermal_control,3);
    const std::array<double,4> reference_reduced={expected[0]/(R*T),expected[1]/(R*T)-expected[0]/(R*T*T),
        expected[2]/(R*T)-2*expected[1]/(R*T*T)+2*expected[0]/(R*T*T*T),
        expected[3]/(R*T)-3*expected[2]/(R*T*T)+6*expected[1]/(R*T*T*T)-6*expected[0]/(R*T*T*T*T)};
    for(std::size_t k=0;k<4;++k) {
        const double scale=std::pow(T,k)/(R*T);
        CHECK(close_with_tolerance(scale*actual.gibbs_derivatives[k],scale*expected[k],1e-10));
        CHECK(close_with_tolerance(scale*targets[k][0],scale*expected[k],1e-10));
        double mapped=0; const double factorial=k<2?1:k==2?2:6;
        for(std::size_t row=0;row<compiled.rows.size();++row) if(compiled.rows[row].kind==RowKind::Chemical)
            mapped+=nu[compiled.rows[row].index]*factorial*reduced.rows[k][row];
        if(k==0) for(std::size_t i=0;i<nu.size();++i) mapped-=nu[i]*std::log(prepared.composition[0][i]/prepared.volume[0]);
        CHECK(close_with_tolerance(std::pow(T,k)*mapped,std::pow(T,k)*reference_reduced[k],1e-10));
        std::printf("thermal-derivative,%zu,actual=%.17g,reference=%.17g,scaled_error=%.17g,reduced=%.17g,reduced_reference=%.17g,reduced_scaled_error=%.17g\n",k,actual.gibbs_derivatives[k],expected[k],scale*(actual.gibbs_derivatives[k]-expected[k]),mapped,reference_reduced[k],std::pow(T,k)*(mapped-reference_reduced[k]));
    }
}
} // namespace equilibrium_checks
