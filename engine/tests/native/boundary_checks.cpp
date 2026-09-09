#include "equilibrium/problem.hpp"
#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium_support.hpp"
#include "recording_support.hpp"

#include "core/numerical_policy.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>

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
    const auto inputs = recording_support::fields({
        recording_support::field("T", 180.0), recording_support::field("x", bubble.feed.values),
        recording_support::field("P_bubble", bubble_result.pressure.front()),
        recording_support::field("P_dew", dew_result.pressure.front()),
        recording_support::field("limit", limit),
    });
    recording_support::compare(
        "pressure-roundtrip", "pressure_ratio", "1",
        dew_result.pressure.front() / bubble_result.pressure.front(), 1.0, inputs,
        "Relationship", "bubble-dew-roundtrip",
        "dew pressure reproduces bubble pressure", "reference_or_one", limit, 0.0,
        std::abs(dew_result.pressure.front() / bubble_result.pressure.front() - 1.0) <= limit,
        "s3-pressure", "bubble-dew", "P", "T,\\mathbf{x},branch", "P_{dew}=P_{bubble}", "branch");
    CHECK(std::abs(dew_result.pressure.front() / bubble_result.pressure.front() - 1.0) <= limit);
    for (std::size_t component = 0; component < 2; ++component) {
        const auto point = std::to_string(component);
        const bool equal = std::abs(
            dew_result.mole_fractions[component] - bubble_result.mole_fractions[component]
        ) <= limit;
        recording_support::compare(
            (std::string("composition-roundtrip/") + point).c_str(), "mole_fraction", "1",
            dew_result.mole_fractions[component], bubble_result.mole_fractions[component],
            recording_support::fields({recording_support::field("T", bubble.temperature.value),
                                       recording_support::field("P_bubble", bubble_result.pressure.front()),
                                       recording_support::field("P_dew", dew_result.pressure.front()),
                                       recording_support::field("feed", bubble.feed.values),
                                       recording_support::field("component", component),
                                       recording_support::field("limit", limit)}),
            "Relationship", "bubble-dew-roundtrip", "dew composition reproduces bubble composition",
            "reference_or_one", limit, 0.0, equal,
            "s3-composition", point.c_str(), "x_i", "T,P_{bubble},P_{dew},feed,component",
            "x_{dew,i}=x_{bubble,i}", "component");
        CHECK(equal);
    }
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
    for (std::size_t step = 1; step < gaps.size(); ++step) {
        const bool contracted = std::isfinite(gaps[step])
            && std::isfinite(gaps[step - 1])
            && gaps[step] <= gaps[step - 1] / 10.0;
        const auto point = std::to_string(step);
        recording_support::compare_custom(
            (std::string("temperature-ladder/") + point).c_str(),
            "temperature_gap", "K", gaps[step], gaps[step - 1],
            recording_support::fields({
                recording_support::field("T", bubble_result.temperature.front()),
                recording_support::field("P", bubble_result.pressure.front()),
                recording_support::field("beta", std::vector<double>{1.0e-2, 1.0e-4, 1.0e-6}),
                recording_support::field("gaps", std::vector<double>{gaps[0], gaps[1], gaps[2]}),
                recording_support::field("step", step),
                recording_support::field("previous_gap", gaps[step - 1]),
            }),
            "Relationship", "boundary-limit", "temperature gap contracts by at least one decade",
            "isfinite(model) && isfinite(reference) && model <= reference / 10", contracted,
            "s4-temperature-ladder", point.c_str(), "\\Delta T", "T,P,\\beta,step,\\Delta T_{previous}",
            "\\Delta T_{k+1}\\leq0.1\\Delta T_k", "step");
        CHECK(contracted);
    }
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
    const bool temperature_matches = std::abs(
        temperature_result.temperature.front() - pressure_result.temperature.front()
    ) <= 1.0e-6;
    recording_support::compare(
        "temperature-roundtrip", "temperature", "K", temperature_result.temperature.front(),
        pressure_result.temperature.front(),
        recording_support::fields({recording_support::field("T", pressure_result.temperature.front()),
                                   recording_support::field("P", pressure_result.pressure.front()),
                                   recording_support::field("x", problem.feed.values)}),
        "Relationship", "pressure-temperature-roundtrip", "temperature solve reproduces pressure solve",
        "reference_or_one", 1.0e-6, 0.0,
        temperature_matches,
        "s5-temperature", "roundtrip", "T", "T,P,\\mathbf{x}", "T_{roundtrip}=T", "state");
    CHECK(temperature_matches);
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
            recording_support::compare(
                (std::string("sensitivity/") + std::to_string(lever_index) + "/"
                 + std::to_string(column)).c_str(),
                "equilibrium_sensitivity", "coordinate", derivative[column], centered,
                recording_support::fields({recording_support::field("T", problem.temperature.value),
                                           recording_support::field("P", problem.pressure.value),
                                           recording_support::field("x", problem.feed.values),
                                           recording_support::field("lever", lever),
                                           recording_support::field("delta", delta),
                                           recording_support::field("coordinate", column)}),
                "Relationship", "finite-difference", "analytic sensitivity matches centered solve",
                "reference_or_one", 0.0, 1.0e-6,
                std::abs(derivative[column] - centered) <= 1.0e-6 * std::max(1.0, std::abs(centered)),
                "s9-sensitivity", (std::to_string(lever_index) + "-" + std::to_string(column)).c_str(),
                "\\partial z/\\partial q", "T,P,\\mathbf{x},\\Delta q", "z_q=\\Delta z/\\Delta q", "lever,coordinate");
            CHECK(std::abs(derivative[column] - centered) <= 1.0e-6
                  * std::max(1.0, std::abs(centered)));
        }
    }

    auto bubble = equilibrium_support::flash(3);
    bubble.phases = {
        equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid),
        equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor,
                                   epcsaft::equilibrium::AmountSpec::pinned(0.0)),
    };
    bubble.temperature.value = 180.0;
    bubble.pressure = epcsaft::equilibrium::FreeValue::variable(1.0e5);
    const auto bubble_result = epcsaft::equilibrium::solve_equilibrium(model, bubble);
    CHECK(bubble_result.success);
    if (!bubble_result.success) return;
    auto pinned = bubble;
    pinned.phases[1].amount = epcsaft::equilibrium::AmountSpec::pinned(1.0e-2);
    pinned.temperature = epcsaft::equilibrium::FreeValue::variable(180.0);
    pinned.pressure = epcsaft::equilibrium::FreeValue::fixed(bubble_result.pressure.front());
    const auto pinned_result = epcsaft::equilibrium::solve_equilibrium(model, pinned, bubble_result);
    CHECK(pinned_result.success);
    if (!pinned_result.success) return;
    const auto pinned_compiled = epcsaft::equilibrium::compile_problem(model, pinned);
    for (std::size_t species = 0; species < 3; ++species) {
        auto higher = pinned;
        auto lower = pinned;
        const double delta = 1.0e-5;
        higher.feed.values[species] += delta;
        lower.feed.values[species] -= delta;
        std::vector<double> lever(5, 0.0);
        lever[2 + species] = 1.0;
        const auto analytic = epcsaft::equilibrium::sensitivity(
            pinned_compiled, pinned_result.coordinates, lever);
        const auto up = epcsaft::equilibrium::solve_equilibrium(model, higher, pinned_result);
        const auto down = epcsaft::equilibrium::solve_equilibrium(model, lower, pinned_result);
        CHECK(up.success && down.success);
        if (!up.success || !down.success) continue;
        for (std::size_t column = 0; column < analytic.size(); ++column) {
            const double centered = (up.coordinates[column] - down.coordinates[column])
                / (2.0 * delta);
            recording_support::compare(
                (std::string("pinned-sensitivity/") + std::to_string(species) + "/"
                 + std::to_string(column)).c_str(),
                "equilibrium_sensitivity", "coordinate", analytic[column], centered,
                recording_support::fields({recording_support::field("T", pinned.temperature.value),
                                           recording_support::field("P", pinned.pressure.value),
                                           recording_support::field("x", pinned.feed.values),
                                           recording_support::field("species", species),
                                           recording_support::field("delta", delta),
                                           recording_support::field("coordinate", column)}),
                "Relationship", "finite-difference", "pinned sensitivity matches centered solve",
                "reference_or_one", 0.0, 1.0e-6,
                std::abs(analytic[column] - centered) <= 1.0e-6 * std::max(1.0, std::abs(centered)),
                "s9-pinned-sensitivity", (std::to_string(species) + "-" + std::to_string(column)).c_str(),
                "\\partial z/\\partial x_i", "T,P,\\mathbf{x},\\Delta x_i", "z_{x_i}=\\Delta z/\\Delta x_i", "species,coordinate");
            CHECK(std::abs(analytic[column] - centered)
                  <= 1.0e-6 * std::max(1.0, std::abs(centered)));
        }
    }

    epcsaft::HelmholtzModel reactive_model(equilibrium_support::ascani_reactive_ternary());
    auto reactive = equilibrium_support::ascani_reactive_problem();
    reactive.phases[1].amount = epcsaft::equilibrium::AmountSpec::pinned(1.0e-2);
    reactive.standard_gibbs_j_per_mol = {1000.0, -1000.0, 0.0};
    reactive.temperature = epcsaft::equilibrium::FreeValue::variable(298.15);
    const auto reactive_free = epcsaft::equilibrium::compile_problem(reactive_model, reactive);
    const auto reactive_start = epcsaft::equilibrium::boundary_start(reactive_free);
    const auto reactive_jacobian = reactive_free.jacobian(reactive_start);
    auto warmer = reactive_start;
    auto colder = reactive_start;
    const std::size_t temperature_column = reactive_free.primal_size - 1;
    const double temperature_step = 1.0e-3;
    warmer[temperature_column] += temperature_step;
    colder[temperature_column] -= temperature_step;
    const auto warm_rows = reactive_free.evaluate(warmer).values;
    const auto cold_rows = reactive_free.evaluate(colder).values;
    for (std::size_t row = 0; row < reactive_free.rows.size(); ++row) {
        const double centered = (warm_rows[row] - cold_rows[row]) / (2.0 * temperature_step);
        recording_support::compare(
            (std::string("reactive-jacobian/") + std::to_string(row)).c_str(),
            "equilibrium_jacobian", "1", reactive_jacobian[row * reactive_start.size() + temperature_column], centered,
            recording_support::fields({recording_support::field("T", 298.15),
                                       recording_support::field("P", 1.0e5),
                                       recording_support::field("x", reactive.feed.values),
                                       recording_support::field("delta_T", temperature_step),
                                       recording_support::field("row", row)}),
            "Relationship", "finite-difference", "reactive Jacobian temperature column matches centered rows",
            "reference_or_one", 0.0, 1.0e-6,
            std::abs(reactive_jacobian[row * reactive_start.size() + temperature_column] - centered)
                <= 1.0e-6 * std::max(1.0, std::abs(centered)),
            "s9-reactive-jacobian", std::to_string(row).c_str(), "J_{r,T}", "T,P,\\mathbf{x},\\Delta T",
            "J_{r,T}=\\Delta r/\\Delta T", "row");
        CHECK(std::abs(reactive_jacobian[row * reactive_start.size() + temperature_column] - centered)
              <= 1.0e-6 * std::max(1.0, std::abs(centered)));
    }
    reactive.temperature = epcsaft::equilibrium::FreeValue::fixed(298.15);
    reactive.pressure = epcsaft::equilibrium::FreeValue::variable();
    const auto reactive_fixed = epcsaft::equilibrium::compile_problem(reactive_model, reactive);
    const auto fixed_start = epcsaft::equilibrium::boundary_start(reactive_fixed);
    const auto reactive_levers = reactive_fixed.lever_jacobian(fixed_start);
    auto warmer_problem = reactive;
    auto colder_problem = reactive;
    warmer_problem.temperature.value += temperature_step;
    colder_problem.temperature.value -= temperature_step;
    const auto warm_system = epcsaft::equilibrium::compile_problem(reactive_model, warmer_problem);
    const auto cold_system = epcsaft::equilibrium::compile_problem(reactive_model, colder_problem);
    const auto fixed_warm_rows = warm_system.evaluate(fixed_start).values;
    const auto fixed_cold_rows = cold_system.evaluate(fixed_start).values;
    const std::size_t lever_width = reactive_levers.size() / reactive_fixed.rows.size();
    for (std::size_t row = 0; row < reactive_fixed.rows.size(); ++row) {
        const double centered = (fixed_warm_rows[row] - fixed_cold_rows[row])
            / (2.0 * temperature_step);
        recording_support::compare(
            (std::string("reactive-lever/") + std::to_string(row)).c_str(), "lever_jacobian", "1",
            reactive_levers[row * lever_width], centered,
            recording_support::fields({recording_support::field("T", 298.15),
                                       recording_support::field("P", 1.0e5),
                                       recording_support::field("x", reactive.feed.values),
                                       recording_support::field("delta_T", temperature_step),
                                       recording_support::field("row", row)}),
            "Relationship", "finite-difference", "lever Jacobian temperature column matches centered rows",
            "reference_or_one", 0.0, 1.0e-6,
            std::abs(reactive_levers[row * lever_width] - centered)
                <= 1.0e-6 * std::max(1.0, std::abs(centered)),
            "s9-reactive-lever", std::to_string(row).c_str(), "J_{r,T}", "T,P,\\mathbf{x},\\Delta T",
            "J_{r,T}=\\Delta r/\\Delta T", "row");
        CHECK(std::abs(reactive_levers[row * lever_width] - centered)
              <= 1.0e-6 * std::max(1.0, std::abs(centered)));
    }

    auto active_parameters = model_support::neutral_ternary();
    const std::vector<epcsaft::ActiveParameter> active = {
        {epcsaft::ParameterFamily::Kij, {1, 3}},
        {epcsaft::ParameterFamily::SegmentCount, {0}},
    };
    epcsaft::HelmholtzModel active_model(active_parameters, active);
    auto active_problem = equilibrium_support::flash(3);
    active_problem.phases = {
        equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor),
        equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid),
    };
    active_problem.temperature.value = 180.0;
    const auto active_result = epcsaft::equilibrium::solve_equilibrium(active_model, active_problem);
    CHECK(active_result.success);
    if (!active_result.success) return;
    const auto active_compiled = epcsaft::equilibrium::compile_problem(active_model, active_problem);
    const auto active_columns = active_compiled.lever_jacobian(active_result.coordinates);
    constexpr std::size_t active_levers = 7;
    for (std::size_t parameter = 0; parameter < 2; ++parameter) {
        std::vector<double> lever(active_levers, 0.0);
        lever[5 + parameter] = 1.0;
        const auto analytic = epcsaft::equilibrium::sensitivity(
            active_compiled, active_result.coordinates, lever);
        auto higher_parameters = active_parameters;
        auto lower_parameters = active_parameters;
        const double theta = parameter == 0 ? active_parameters.k_ij[1]
            : active_parameters.segment_count[0];
        const double delta = 1.0e-4 * std::max(std::abs(theta), 1.0);
        if (parameter == 0) {
            higher_parameters.k_ij[1] += delta;
            higher_parameters.k_ij[3] += delta;
            lower_parameters.k_ij[1] -= delta;
            lower_parameters.k_ij[3] -= delta;
        } else {
            higher_parameters.segment_count[0] += delta;
            lower_parameters.segment_count[0] -= delta;
        }
        const auto up = epcsaft::equilibrium::solve_equilibrium(
            epcsaft::HelmholtzModel(higher_parameters, active), active_problem, active_result);
        const auto down = epcsaft::equilibrium::solve_equilibrium(
            epcsaft::HelmholtzModel(lower_parameters, active), active_problem, active_result);
        CHECK(up.success && down.success);
        for (std::size_t column = 0; up.success && down.success && column < analytic.size(); ++column) {
            const double centered = (up.coordinates[column] - down.coordinates[column])
                / (2.0 * delta);
            recording_support::compare(
                (std::string("parameter-sensitivity/") + std::to_string(parameter) + "/"
                 + std::to_string(column)).c_str(), "equilibrium_sensitivity", "coordinate",
                analytic[column], centered,
                recording_support::fields({recording_support::field("T", active_problem.temperature.value),
                                           recording_support::field("P", active_problem.pressure.value),
                                           recording_support::field("x", active_problem.feed.values),
                                           recording_support::field("parameter", parameter),
                                           recording_support::field("theta", theta),
                                           recording_support::field("delta", delta),
                                           recording_support::field("coordinate", column)}),
                "Relationship", "finite-difference", "active parameter sensitivity matches centered solve",
                "reference_or_one", 0.0, 1.0e-6,
                std::abs(analytic[column] - centered) <= 1.0e-6 * std::max(1.0, std::abs(centered)),
                "s9-parameter-sensitivity", (std::to_string(parameter) + "-" + std::to_string(column)).c_str(),
                "\\partial z/\\partial\\theta", "T,P,\\mathbf{x},\\theta,\\Delta\\theta",
                "z_\\theta=\\Delta z/\\Delta\\theta", "parameter,coordinate");
            CHECK(std::abs(analytic[column] - centered)
                  <= 1.0e-6 * std::max(1.0, std::abs(centered)));
        }
        for (std::size_t row = 0; row < active_compiled.rows.size(); ++row)
            if (active_compiled.rows[row].kind == epcsaft::equilibrium::RowKind::Moiety)
                CHECK(active_columns[active_levers * row + 5 + parameter] == 0.0);
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
        for (std::size_t row = 0; row < values.size(); ++row) {
            const double centered = (up[row] - down[row]) / (2.0 * delta);
            const auto point_id = std::to_string(row) + "-" + std::to_string(column);
            recording_support::compare(
                (std::string("jacobian/") + point_id).c_str(), "equilibrium_jacobian", "1",
                jacobian[row * point.size() + column], centered,
                recording_support::fields({recording_support::field("T", 300.0),
                                           recording_support::field("P", 1.0e5),
                                           recording_support::field("x", std::vector<double>{0.5, 0.5}),
                                           recording_support::field("row", row),
                                           recording_support::field("column", column),
                                           recording_support::field("delta", delta)}),
                "Relationship", "finite-difference", "equilibrium Jacobian matches centered residual derivative",
                "reference_or_one", 0.0, 1.0e-6,
                std::abs(jacobian[row * point.size() + column] - centered)
                    <= 1.0e-6 * std::max(1.0, std::abs(centered)),
                "s10-jacobian", point_id.c_str(), "J", "T,P,\\mathbf{x},\\Delta z",
                "J_{ij}=\\Delta r_i/\\Delta z_j", "row,col");
            CHECK(std::abs(jacobian[row * point.size() + column] - centered)
                  <= 1.0e-6 * std::max(1.0, std::abs(centered)));
        }
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
    recording_support::compare(
        "pinned-pressure", "pressure_ratio", "1",
        vapor_pinned.pressure.front() / liquid_pinned.pressure.front(), 1.0,
        recording_support::fields({recording_support::field("T", 250.0),
                                   recording_support::field("P", 1.3e6),
                                   recording_support::field("limit", limit),
                                   recording_support::field("x", std::vector<double>{1.0})}),
        "Relationship", "pinned-phase-switch", "switching pinned phase preserves pressure",
        "reference_or_one", limit, 0.0,
        std::abs(vapor_pinned.pressure.front() / liquid_pinned.pressure.front() - 1.0) <= limit,
        "s14-pressure", "phase-switch", "P_{vapor}/P_{liquid}", "T,P,\\mathbf{x},branch",
        "P_{vapor}=P_{liquid}", "branch");
    CHECK(std::abs(vapor_pinned.pressure.front() / liquid_pinned.pressure.front() - 1.0) <= limit);
}

void check_s15c() {
    epcsaft::HelmholtzModel model(model_support::neutral_binary());
    auto problem = equilibrium_support::boundary();
    problem.feed.values = {0.5, 0.5};
    problem.temperature.value = 280.0;
    const auto result = epcsaft::equilibrium::solve_equilibrium(model, problem);
    recording_support::compare(
        "coincident-failure", "success", "bool", result.success ? 1.0 : 0.0, 0.0,
        recording_support::fields({recording_support::field("T", 280.0),
                                   recording_support::field("P", 1.0e5),
                                   recording_support::field("x", problem.feed.values),
                                   recording_support::field("topology", "boundary")}),
        "Relationship", "boundary-admission", "coincident phases are rejected",
        "reference_or_one", 0.0, 0.0, !result.success, "s15c-failure", "boundary",
        "success", "T,P,\\mathbf{x},topology", "success=0 for coincident phases", "topology");
    CHECK(!result.success);
    CHECK(!result.has_state());
    CHECK(result.message.find("phases_coincide") != std::string::npos);
    CHECK(result.message.find("distance") != std::string::npos);
}

}  // namespace

int main(int argc, char** argv) {
    return recording_support::run(argc, argv, {
        {{"native.boundary.s3", "Boundary_S3", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "pressure-roundtrip,composition-roundtrip/*", "engine/tests/native/boundary_checks.cpp", "boundary_checks", "S3"}, check_s3},
        {{"native.boundary.s4", "Boundary_S4", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "temperature-ladder/*", "engine/tests/native/boundary_checks.cpp", "boundary_checks", "S4", "temperature-ladder/*:2"}, check_s4},
        {{"native.boundary.s5", "Boundary_S5", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "temperature-roundtrip", "engine/tests/native/boundary_checks.cpp", "boundary_checks", "S5"}, check_s5},
        {{"native.boundary.s9", "Boundary_S9", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "sensitivity/*,pinned-sensitivity/*,reactive-jacobian/*,reactive-lever/*,parameter-sensitivity/*", "engine/tests/native/boundary_checks.cpp", "boundary_checks", "S9"}, check_s9},
        {{"native.boundary.s10", "Boundary_S10", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "jacobian/*", "engine/tests/native/boundary_checks.cpp", "boundary_checks", "S10"}, check_s10},
        {{"native.boundary.s14", "Boundary_S14", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "pinned-pressure", "engine/tests/native/boundary_checks.cpp", "boundary_checks", "S14"}, check_s14},
        {{"native.boundary.s15c", "Boundary_S15c", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "coincident-failure", "engine/tests/native/boundary_checks.cpp", "boundary_checks", "S15c"}, check_s15c},
    });
}
