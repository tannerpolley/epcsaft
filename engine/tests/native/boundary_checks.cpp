#include "equilibrium/problem.hpp"
#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium_support.hpp"
#include "recording_support.hpp"

#include "core/numerical_policy.hpp"
#include "equilibrium/stoichiometry.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <optional>
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
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    auto bubble = equilibrium_support::boundary();
    bubble.feed.values = {0.5, 0.5};
    bubble.temperature.value = 180.0;
    bubble.pressure.value = 1.0e5;
    const auto bubble_result = epcsaft::equilibrium::solve_equilibrium(model, bubble);
    CHECK(bubble_result.success);
    auto dew = bubble;
    dew.feed.values = {bubble_result.mole_fractions[2], bubble_result.mole_fractions[3]};
    dew.pressure.value = bubble_result.pressure.front();
    dew.phases[0].amount = epcsaft::equilibrium::AmountSpec::pinned(0.0);
    dew.phases[1].amount = epcsaft::equilibrium::AmountSpec::free();
    const auto dew_result = epcsaft::equilibrium::solve_equilibrium(model, dew);
    CHECK(dew_result.success);
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
    }
}

void check_s4() {
    epcsaft::HelmholtzModel model(properties_support::neutral_ternary());
    auto bubble = equilibrium_support::flash(3);
    bubble.phases = {equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid),
                     equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor,
                                                epcsaft::equilibrium::AmountSpec::pinned(0.0))};
    bubble.temperature.value = 180.0;
    bubble.pressure = epcsaft::equilibrium::FreeValue::variable(1.0e5);
    const auto bubble_result = epcsaft::equilibrium::solve_equilibrium(model, bubble);
    CHECK(bubble_result.success);
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
    }
}

void check_s5() {
    // A TP start must satisfy the feed balance before the general driver iterates.
    epcsaft::HelmholtzModel binary(properties_support::neutral_binary());
    auto flash = equilibrium_support::flash(2);
    flash.phases = {equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor),
                    equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid)};
    flash.temperature.value = 180.0;
    flash.pressure.value = 1.0e6;
    flash.feed.values = {5.0, 5.0};
    const auto compiled = epcsaft::equilibrium::compile_problem(binary, flash);
    const auto seeded = compiled.prepare(epcsaft::equilibrium::boundary_start(compiled));
    CHECK(seeded.amount[0] > 0.0 && seeded.amount[1] > 0.0);
    for (std::size_t i = 0; i < 2; ++i) {
        CHECK(std::abs(seeded.amount[0] * seeded.composition[0][i]
                     + seeded.amount[1] * seeded.composition[1][i] - 5.0) < 1.0e-12);
        CHECK(std::abs(seeded.composition[i][0] + seeded.composition[i][1] - 1.0) < 1.0e-12);
    }
    // A nonbracketed K estimate retains the original equal-amount fallback.
    flash.feed.values = {9.99, 0.01};
    const auto edge = epcsaft::equilibrium::compile_problem(binary, flash);
    const auto fallback = edge.prepare(epcsaft::equilibrium::boundary_start(edge));
    CHECK(fallback.amount[0] == 5.0 && fallback.amount[1] == 5.0);
    epcsaft::HelmholtzModel model(properties_support::neutral_ternary());
    auto problem = equilibrium_support::flash(3);
    problem.phases = {equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid),
                      equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor,
                                                 epcsaft::equilibrium::AmountSpec::pinned(0.0))};
    problem.temperature.value = 180.0;
    problem.pressure = epcsaft::equilibrium::FreeValue::variable(1.0e5);
    const auto pressure_result = epcsaft::equilibrium::solve_equilibrium(model, problem);
    CHECK(pressure_result.success);
    problem.temperature = epcsaft::equilibrium::FreeValue::variable(
        pressure_result.temperature.front());
    problem.pressure = epcsaft::equilibrium::FreeValue::fixed(
        pressure_result.pressure.front());
    const auto temperature_result = epcsaft::equilibrium::solve_equilibrium(
        model, problem, pressure_result);
    CHECK(temperature_result.success);
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
}

void check_s9() {
    epcsaft::HelmholtzModel model(properties_support::neutral_ternary());
    auto problem = equilibrium_support::flash(3);
    problem.phases = {equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor),
                      equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid)};
    problem.temperature.value = 180.0;
    const auto result = epcsaft::equilibrium::solve_equilibrium(model, problem);
    CHECK(result.success);
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
    auto pinned = bubble;
    pinned.phases[1].amount = epcsaft::equilibrium::AmountSpec::pinned(1.0e-2);
    pinned.temperature = epcsaft::equilibrium::FreeValue::variable(180.0);
    pinned.pressure = epcsaft::equilibrium::FreeValue::fixed(bubble_result.pressure.front());
    const auto pinned_result = epcsaft::equilibrium::solve_equilibrium(model, pinned, bubble_result);
    CHECK(pinned_result.success);
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
        }
    }

    epcsaft::HelmholtzModel reactive_model(equilibrium_support::synthetic_reactive_ternary());
    auto reactive = equilibrium_support::synthetic_reactive_problem();
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
    }

    auto active_parameters = properties_support::neutral_ternary();
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
        }
        for (std::size_t row = 0; row < active_compiled.rows.size(); ++row)
            if (active_compiled.rows[row].kind == epcsaft::equilibrium::RowKind::Moiety)
                CHECK(active_columns[active_levers * row + 5 + parameter] == 0.0);
    }
}

void check_ionic_packing_derivative() {
    const auto parameters = properties_support::aqueous_electrolyte();
    epcsaft::HelmholtzModel model(parameters, {{epcsaft::ParameterFamily::PackingDiameter, {1}}});
    auto problem = equilibrium_support::flash(3);
    problem.feed.values = {0.98, 0.01, 0.01};
    const auto central = epcsaft::equilibrium::solve_equilibrium(model, problem);
    CHECK(central.success);
    const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
    const auto columns = compiled.lever_jacobian(central.coordinates);
    const std::size_t width = 2 + parameters.component_count() + model.active_parameters().size();
    // The two-step study's largest absolute stencil discrepancy was 1.35e-9 / Angstrom.
    // A 2e-9 absolute criterion covers that rounding/refinement error, not the missing volume chain.
    constexpr double tolerance = 2.0e-9;
    for (const double step : {1.0e-4, 1.0e-5}) {
        auto higher = parameters, lower = parameters;
        higher.packing_diameter_angstrom[1] += step;
        lower.packing_diameter_angstrom[1] -= step;
        epcsaft::HelmholtzModel up_model(higher), down_model(lower);
        const auto up = epcsaft::equilibrium::compile_problem(up_model, problem).evaluate(central.coordinates).values;
        const auto down = epcsaft::equilibrium::compile_problem(down_model, problem).evaluate(central.coordinates).values;
        for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
            const auto kind = compiled.rows[row].kind;
            if (kind != epcsaft::equilibrium::RowKind::Chemical && kind != epcsaft::equilibrium::RowKind::Pressure) continue;
            const double exact = columns[row * width + width - 1];
            const double centered = (up[row] - down[row]) / (2.0 * step);
            const bool agrees = std::abs(exact - centered) <= tolerance;
            const auto id = "ionic-packing-derivative/" + recording_support::number(step) + "/" + std::to_string(row);
            recording_support::compare(id.c_str(), "lever_jacobian", "1/angstrom", exact, centered,
                recording_support::fields({recording_support::field("T", 300.0), recording_support::field("P", 1.0e5),
                    recording_support::field("feed", problem.feed.values), recording_support::field("step", step),
                    recording_support::field("row", row)}), "Relationship", "centered-parameter-stencil",
                "ionic packing-diameter column includes the volume chain", "absolute", tolerance, 0.0, agrees);
        }
    }
}

void check_s10() {
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
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
        }
    }
}

void check_s14() {
    epcsaft::HelmholtzModel model(properties_support::ethane());
    auto problem = equilibrium_support::flash(1);
    problem.phases = {equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid),
                      equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor,
                                                 epcsaft::equilibrium::AmountSpec::pinned(0.0))};
    problem.temperature.value = 250.0;
    problem.pressure = epcsaft::equilibrium::FreeValue::variable(1.3e6);
    const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
    CHECK(compiled.rank == 1 && compiled.generic_nullity == 0);
    CHECK(compiled.charge_row_count == 0);
    for (const auto& row : compiled.rows) CHECK(row.kind != epcsaft::equilibrium::RowKind::Electro);
    for (const auto& row : compiled.rows) std::printf("S14 row %s\n", epcsaft::equilibrium::row_kind_name(row.kind));
    const auto vapor_pinned = epcsaft::equilibrium::solve_equilibrium(model, problem);
    CHECK(vapor_pinned.success);
    problem.phases[0].amount = epcsaft::equilibrium::AmountSpec::pinned(0.0);
    problem.phases[1].amount = epcsaft::equilibrium::AmountSpec::free();
    const auto liquid_pinned = epcsaft::equilibrium::solve_equilibrium(model, problem, vapor_pinned);
    CHECK(liquid_pinned.success);
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
}

void check_s15c() {
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
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
    CHECK(!result.has_state());
    CHECK(result.message.find("phases_coincide") != std::string::npos);
    CHECK(result.message.find("distance") != std::string::npos);
}

// ------------------------------------------------------------------ E4 continuation

namespace e4 {

using epcsaft::equilibrium::AmountSpec;
using epcsaft::equilibrium::ContinuationOptions;
using epcsaft::equilibrium::ContinuationRequest;
using epcsaft::equilibrium::ContinuationResult;
using epcsaft::equilibrium::ContinuationStop;
using epcsaft::equilibrium::EquilibriumProblem;
using epcsaft::equilibrium::EquilibriumResult;
using epcsaft::equilibrium::FreeValue;
using epcsaft::equilibrium::PhaseKind;

ContinuationOptions options(double initial, double minimum, std::size_t steps, std::size_t retries) {
    ContinuationOptions result;
    result.initial_step = initial; result.minimum_step = minimum;
    result.max_steps = steps; result.max_retries_per_step = retries;
    return result;
}

std::string summary(const ContinuationResult& result) {
    return recording_support::fields({
        recording_support::field("stop", epcsaft::equilibrium::continuation_stop_name(result.stop)),
        recording_support::field("path_values", result.path_values),
        recording_support::field("reached_path", result.reached_path),
        recording_support::field("attempted_path", epcsaft::core::optional_number(result.attempted_path)),
        recording_support::field("attempts", result.attempts),
        recording_support::field("accepted_steps", result.accepted_steps),
        recording_support::field("retries", result.retries),
        recording_support::field("event_phase", result.event_phase.has_value() ? std::to_string(*result.event_phase) : "absent"),
        recording_support::field("cause", result.cause_stop.has_value() ? epcsaft::equilibrium::continuation_stop_name(*result.cause_stop) : "absent"),
        recording_support::field("message", result.message)});
}

bool record_bool(const char* id, const char* group, const char* claim, bool verdict, const std::string& inputs) {
    recording_support::compare(id, "verdict", "bool", verdict ? 1.0 : 0.0, 1.0, inputs, "Relationship",
        "continuation", claim, "absolute", 0.0, 0.0, verdict, group, id, "verdict", "inputs", claim, "case");
    return verdict;
}

// Corrected-phase branch identity recomputed through the public density owner: the corrected
// density and the anchored locally unique root share the anchor window; slope positive.
void branch_rows(const char* group, const char* label, const epcsaft::HelmholtzModel& model,
                 const EquilibriumProblem& problem, const ContinuationResult& result) {
    const std::size_t components = model.parameters().component_count();
    std::vector<double> anchors;
    for (std::size_t k = 0; k < result.states.size(); ++k) {
        const auto& state = result.states[k];
        for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
            const std::vector<double> x(state.mole_fractions.begin() + static_cast<std::ptrdiff_t>(phase * components),
                                        state.mole_fractions.begin() + static_cast<std::ptrdiff_t>((phase + 1) * components));
            const auto branch = problem.phases[phase].kind == PhaseKind::Vapor ? epcsaft::DensityBranch::Vapor : epcsaft::DensityBranch::Liquid;
            const double corrected = state.molar_densities[phase];
            const double anchor = k == 0
                ? epcsaft::solve_density(model, state.temperature[0], state.pressure[0], x, branch, std::nullopt).molar_density
                : anchors[phase];
            const auto root = epcsaft::solve_density(model, state.temperature[0], state.pressure[0], x, branch, anchor);
            const auto physical = model.state(state.temperature[0], corrected, x);
            const double displacement = std::abs(std::log(corrected / anchor));
            const double selected_displacement = std::abs(std::log(root.molar_density / anchor));
            const double pressure_row = [&] {
                double worst = 0.0;
                for (std::size_t row = 0; row < state.rows.size(); ++row)
                    if (state.rows[row].kind == epcsaft::equilibrium::RowKind::Pressure && state.rows[row].phase == phase)
                        worst = std::max(worst, std::abs(state.classified_residuals[row]));
                return worst;
            }();
            const bool verdict = std::isfinite(corrected) && corrected > 0.0 && physical.pressure_density_derivative > 0.0
                && root.convergence.residual <= epcsaft::policy::DENSITY_ROOT_RESIDUAL.value
                && displacement <= epcsaft::policy::DENSITY_ANCHOR_LOG_WINDOW.value
                && selected_displacement <= epcsaft::policy::DENSITY_ANCHOR_LOG_WINDOW.value
                && pressure_row <= problem.row_tolerance && state.requested_tolerance_met.value_or(false);
            const auto id = std::string(group) + "/branch/" + label + "/" + std::to_string(k) + "/" + std::to_string(phase);
            recording_support::compare_custom(id.c_str(), "corrected_displacement", "1", displacement,
                epcsaft::policy::DENSITY_ANCHOR_LOG_WINDOW.value,
                recording_support::fields({recording_support::field("q", result.path_values[k]),
                    recording_support::field("phase", phase), recording_support::field("selected_density", root.molar_density),
                    recording_support::field("corrected_density", corrected), recording_support::field("anchor", anchor),
                    recording_support::field("selected_displacement", selected_displacement),
                    recording_support::field("corrected_slope", physical.pressure_density_derivative),
                    recording_support::field("branch", root.branch == epcsaft::DensityBranch::Vapor ? "vapor" : root.branch == epcsaft::DensityBranch::Liquid ? "liquid" : "unspecified"),
                    recording_support::field("stable_root_count", root.stable_root_count.has_value() ? std::to_string(*root.stable_root_count) : "absent"),
                    recording_support::field("density_residual", root.convergence.residual),
                    recording_support::field("pressure_row", pressure_row),
                    recording_support::field("row_residual", state.convergence.residual)}),
                "Relationship", "anchored-branch", "corrected and selected densities share the anchored locally unique window",
                "model <= reference && selected_displacement <= reference && slope > 0 && pressure_row <= row_tolerance", verdict, group, id.c_str(), "|\\\\ln(\\\\rho/\\\\rho_a)|", "q,phase", "|\\\\ln(\\\\rho/\\\\rho_a)|\\\\leq w", "q,phase");
        }
        anchors = state.molar_densities;
    }
}

std::size_t attempt_rows(const char* group, const epcsaft::core::Trace& trace) {
    std::size_t rows = 0;
    for (std::size_t span = 0; span < trace.spans.size(); ++span) {
        const std::string name = trace.spans[span].name;
        if (name != "equilibrium_continuation_attempt" && name != "equilibrium_continuation_initial") continue;
        double q = std::numeric_limits<double>::quiet_NaN(), accepted = 0.0, present = 0.0, projected = 0.0;
        std::string cause = "absent";
        std::optional<double> residual;
        for (const auto& event : trace.events) {
            if (event.span != span) continue;
            const std::string label = event.name;
            if (label == "path_value") q = event.values.front();
            else if (label == "accepted") accepted = event.values.front();
            else if (label == "driver_result_present") present = event.values.front();
            else if (label == "row_residual") residual = event.values.front();
            else if (label == "predictor_projected") projected = event.values.front();
            else if (label != "selected_density" && label != "corrected_density" && label != "anchor_density"
                     && label != "corrected_displacement" && label != "corrected_slope" && label != "density_residual") cause = label;
        }
        const auto id = std::string(group) + "/attempt/" + std::to_string(rows);
        recording_support::observe(id.c_str(), "path_value", "1", q, std::nullopt,
            recording_support::fields({recording_support::field("span", name), recording_support::field("accepted", accepted),
                recording_support::field("driver_result_present", present), recording_support::field("cause", cause),
                recording_support::field("predictor_projected", projected),
                recording_support::field("row_residual", residual.has_value() ? recording_support::number(*residual) : "unmeasured")}),
            "Relationship", "continuation-attempt", "retained attempt from the existing Core trace", group, id.c_str());
        ++rows;
    }
    return rows;
}

}  // namespace e4

void check_e4_a() {
    using namespace e4;
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    ContinuationRequest request;
    request.start_problem = equilibrium_support::bubble_at(180.0);
    request.target_problem = equilibrium_support::bubble_at(181.0);
    request.options = options(0.25, 0.25, 4, 2);
    epcsaft::core::Capture capture;
    const auto result = epcsaft::equilibrium::continue_equilibrium(model, request);
    const std::vector<double> expected = {0.0, 0.25, 0.5, 0.75, 1.0};
    bool path_ok = result.success && result.path_values.size() == 5 && result.states.size() == 5;
    for (std::size_t k = 0; path_ok && k < 5; ++k) path_ok = std::abs(result.path_values[k] - expected[k]) <= 1.0e-15;
    record_bool("e4-a/path", "e4-a", "declared bubble branch reaches 181 K at the declared path values",
  path_ok && !result.event_phase.has_value() && !result.event_pair.has_value(), summary(result));
    for (std::size_t k = 0; k < 5; ++k) {
        const auto cold = epcsaft::equilibrium::solve_equilibrium(model, equilibrium_support::bubble_at(180.0 + expected[k]));
        const auto& state = result.states[k];
        const double pressure_error = std::abs(state.pressure.front() / cold.pressure.front() - 1.0);
        double composition_error = 0.0;
        for (std::size_t i = 0; i < state.mole_fractions.size(); ++i)
            composition_error = std::max(composition_error, std::abs(state.mole_fractions[i] - cold.mole_fractions[i]));
        const bool verdict = cold.success && pressure_error <= 1.0e-8 && composition_error <= 1.0e-8;
        const auto id = std::string("e4-a/cold/") + std::to_string(k);
        recording_support::compare(id.c_str(), "pressure_relative_error", "1", pressure_error, 0.0,
            recording_support::fields({recording_support::field("q", expected[k]), recording_support::field("T", 180.0 + expected[k]),
                recording_support::field("P_path", state.pressure.front()), recording_support::field("P_cold", cold.pressure.front()),
                recording_support::field("composition_error", composition_error), recording_support::field("x_path", state.mole_fractions)}),
            "Relationship", "branch-identity", "continuation state equals the independent cold solve",
            "strict_absolute", 1.0e-8, 0.0, verdict, "e4-a", id.c_str(), "\\\\Delta P/P", "q,T", "\\\\Delta P/P\\\\leq10^{-8}", "q");
    }
    branch_rows("e4-a", "bubble", model, request.start_problem, result);
    attempt_rows("e4-a", capture.trace());
}

void check_e4_b() {
    using namespace e4;
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    const auto start = equilibrium_support::interior(model);
    CHECK(start.start.success && start.start.requested_tolerance_met.value_or(false));
    const std::vector<double> x_l = {start.start.mole_fractions[0], start.start.mole_fractions[1]};
    const std::vector<double> x_v = {start.start.mole_fractions[2], start.start.mole_fractions[3]};
    ContinuationRequest request;
    request.start_problem = start.problem;
    request.target_problem = start.problem;
    request.target_problem.feed.values = {0.4 * x_l[0] + 0.6 * x_v[0], 0.4 * x_l[1] + 0.6 * x_v[1]};
    request.start_result = start.start;
    request.options = options(0.25, 0.25, 4, 2);
    epcsaft::core::Capture capture;
    const auto result = epcsaft::equilibrium::continue_equilibrium(model, request);
    record_bool("e4-b/path", "e4-b", "interior tie line follows the feed lever to the target", result.success, summary(result));
    for (std::size_t k = 0; k < result.states.size(); ++k) {
        const double q = result.path_values[k];
        const auto& state = result.states[k];
        std::vector<double> z_input(2), z_lever(2, 0.0);
        for (std::size_t i = 0; i < 2; ++i) z_input[i] = (1.0 - q) * start.feed[i] + q * request.target_problem.feed.values[i];
        double amount_total = 0.0;
        for (std::size_t phase = 0; phase < 2; ++phase) {
            amount_total += state.amounts[phase];
            for (std::size_t i = 0; i < 2; ++i) z_lever[i] += state.amounts[phase] * state.mole_fractions[phase * 2 + i];
        }
        double lever_error = 0.0, composition_error = 0.0, amount_error = 0.0;
        for (std::size_t i = 0; i < 2; ++i) {
            lever_error = std::max(lever_error, std::abs(z_lever[i] / amount_total - z_input[i]));
            composition_error = std::max({composition_error, std::abs(state.mole_fractions[i] - x_l[i]), std::abs(state.mole_fractions[2 + i] - x_v[i])});
        }
        const double feed_total = z_input[0] + z_input[1];
        const double total_error = std::abs(amount_total - feed_total) / std::max(1.0, feed_total);
        amount_error = std::max(std::abs(state.amounts[0] - ((1.0 - q) * start.start.amounts[0] + 0.4 * q)),
                                std::abs(state.amounts[1] - ((1.0 - q) * start.start.amounts[1] + 0.6 * q)));
        const bool verdict = lever_error <= 1.0e-8 && total_error <= 1.0e-10 && composition_error <= 1.0e-8 && amount_error <= 1.0e-8;
        const auto id = std::string("e4-b/lever/") + std::to_string(k);
        recording_support::compare(id.c_str(), "lever_rule_error", "1", lever_error, 0.0,
            recording_support::fields({recording_support::field("q", q), recording_support::field("z_input", z_input),
                recording_support::field("amounts", state.amounts), recording_support::field("x", state.mole_fractions),
                recording_support::field("total_error", total_error), recording_support::field("composition_error", composition_error),
                recording_support::field("amount_error", amount_error)}),
            "Relationship", "lever-rule", "returned amounts and compositions reconstruct the input feed",
            "strict_absolute", 1.0e-8, 0.0, verdict, "e4-b", id.c_str(), "\\\\Delta z", "q", "\\\\Delta z\\\\leq10^{-8}", "q");
    }
    branch_rows("e4-b", "tie-line", model, start.problem, result);
    attempt_rows("e4-b", capture.trace());
}

void check_e4_c() {
    using namespace e4;
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    const auto start = equilibrium_support::interior(model);
    ContinuationRequest request;
    request.start_problem = start.problem;
    request.target_problem = start.problem;
    request.target_problem.pressure = FreeValue::fixed(2.0 * start.bubble_pressure);
    request.start_result = start.start;
    request.options = options(0.002, 0.000125, 8192, 4);
    epcsaft::core::Capture capture;
    const auto result = epcsaft::equilibrium::continue_equilibrium(model, request);
    bool distinct = true;
    for (std::size_t k = 1; k < result.states.size(); ++k)
        distinct = distinct && result.states[k].amounts != result.states[k - 1].amounts;
    record_bool("e4-c/event", "e4-c", "the vapor phase disappears before twice the bubble pressure",
  !result.success && result.stop == ContinuationStop::phase_disappeared && result.event_phase == std::optional<std::size_t>(1)
      && result.attempted_path > result.reached_path && distinct && result.states.size() == result.path_values.size(),
  recording_support::fields({summary(result), recording_support::field("P_bubble", start.bubble_pressure),
      recording_support::field("feed", start.feed), recording_support::field("start_coordinates", start.start.coordinates),
      recording_support::field("states", result.states.size())}));
    branch_rows("e4-c", "prefix", model, start.problem, result);
    attempt_rows("e4-c", capture.trace());
}

void check_e4_d() {
    using namespace e4;
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    ContinuationRequest request;
    request.start_problem = equilibrium_support::bubble_at(180.0);
    request.target_problem = equilibrium_support::bubble_at(280.0);
    request.options = options(0.001, 0.0000625, 20000, 4);
    epcsaft::core::Capture capture;
    const auto result = epcsaft::equilibrium::continue_equilibrium(model, request);
    // The bubble branch passes through the mixture critical point: the pair's identities swap
    // between consecutive accepted states, which the controller reports as the pair's coalescence
    // (owner decision 2026-09-16). A step or window limit is a failed event row here.
    const bool coalesced = result.stop == ContinuationStop::phases_coalesced && result.event_pair == std::make_pair(std::size_t{0}, std::size_t{1});
    record_bool("e4-d/event", "e4-d", "the bubble branch ends at the named coalescence of L and V, never at the target",
  !result.success && coalesced && !result.states.empty() && result.reached_path < 1.0 && result.last_result.has_value(),
  recording_support::fields({summary(result), recording_support::field("last_T", 180.0 + 100.0 * result.reached_path),
      recording_support::field("states", result.states.size())}));
    branch_rows("e4-d", "bubble", model, request.start_problem, result);
    // The swap itself: every nonzero component of the inter-phase difference vector (delta ln eta,
    // delta x_i) changes sign between the last accepted state and the rejected swapped state.
    if (!result.states.empty() && result.last_result.has_value() && result.last_result->has_state()) {
        const auto& before = result.states.back();
        const auto& after = *result.last_result;
        std::vector<double> before_delta{std::log(before.packing_fractions[0] / before.packing_fractions[1])};
        std::vector<double> after_delta{std::log(after.packing_fractions[0] / after.packing_fractions[1])};
        for (std::size_t i = 0; i < 2; ++i) {
            before_delta.push_back(before.mole_fractions[i] - before.mole_fractions[2 + i]);
            after_delta.push_back(after.mole_fractions[i] - after.mole_fractions[2 + i]);
        }
        double worst = -std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < before_delta.size(); ++i) worst = std::max(worst, before_delta[i] * after_delta[i]);
        recording_support::compare_custom("e4-d/swap", "sign_product_max", "1", worst, 0.0,
            recording_support::fields({recording_support::field("T_before", before.temperature.front()),
                recording_support::field("T_after", after.temperature.front()),
                recording_support::field("before_delta", before_delta), recording_support::field("after_delta", after_delta),
                recording_support::field("rho_before", before.molar_densities), recording_support::field("rho_after", after.molar_densities),
                recording_support::field("x_before", before.mole_fractions), recording_support::field("x_after", after.mole_fractions),
                recording_support::field("row_residual_after", after.convergence.residual)}),
            "Relationship", "phase-swap", "every nonzero inter-phase difference component changes sign across the critical point",
            "model < reference", worst < 0.0, "e4-d", "swap", "\\max_i \\Delta_i^- \\Delta_i^+", "component", "\\Delta_i^- \\Delta_i^+ < 0", "component");
    }
    // Retained degeneration measurements at the last accepted state before the swap: phase
    // distance, corrected slopes and densities.
    if (!result.states.empty()) {
        const auto& last = result.states.back();
        const auto compiled = epcsaft::equilibrium::compile_problem(model, equilibrium_support::bubble_at(180.0 + 100.0 * result.reached_path));
        const auto topology = epcsaft::equilibrium::check_topology(compiled, last.coordinates);
        std::vector<double> slopes;
        for (std::size_t phase = 0; phase < 2; ++phase) {
            const std::vector<double> x(last.mole_fractions.begin() + static_cast<std::ptrdiff_t>(2 * phase),
                                        last.mole_fractions.begin() + static_cast<std::ptrdiff_t>(2 * phase + 2));
            slopes.push_back(model.state(last.temperature.front(), last.molar_densities[phase], x).pressure_density_derivative);
        }
        recording_support::observe("e4-d/last-state", "phase_distance", "1", topology.distances.empty() ? 0.0 : topology.distances.front(), std::nullopt,
            recording_support::fields({recording_support::field("T", last.temperature.front()), recording_support::field("P", last.pressure.front()),
                recording_support::field("x", last.mole_fractions), recording_support::field("rho", last.molar_densities),
                recording_support::field("corrected_slopes", slopes), recording_support::field("row_residual", last.convergence.residual),
                recording_support::field("coincidence_threshold", epcsaft::policy::PHASE_COINCIDENCE_DISTANCE.value)}),
            "Relationship", "anchored-branch", "degeneration measurements at the last accepted state before the swap",
            "e4-d", "last-state");
    }
    attempt_rows("e4-d", capture.trace());
}

void check_e4_e() {
    using namespace e4;
    epcsaft::HelmholtzModel model(equilibrium_support::synthetic_reactive_ternary());
    auto problem = equilibrium_support::synthetic_reactive_problem();
    const auto first = epcsaft::equilibrium::solve_equilibrium(model, problem);
    recording_support::observe("e4-e/start", "row_residual", "1", first.convergence.residual, std::nullopt,
        recording_support::fields({recording_support::field("solver_status", first.solver_status),
            recording_support::field("stop", epcsaft::core::stop_name(first.convergence.stop)),
            recording_support::field("floor", epcsaft::core::optional_number(first.convergence.floor)),
            recording_support::field("row_tolerance", problem.row_tolerance),
            recording_support::field("requested_tolerance_met", first.requested_tolerance_met.has_value() ? (*first.requested_tolerance_met ? "true" : "false") : "absent"),
            recording_support::field("message", first.message)}),
        "Relationship", "reactive-continuation", "accepted reactive two-liquid start at the caller criterion", "e4-e", "start");
    CHECK(first.success && first.requested_tolerance_met.value_or(false));
    std::vector<double> x1(first.mole_fractions.begin(), first.mole_fractions.begin() + 3);
    std::vector<double> x2(first.mole_fractions.begin() + 3, first.mole_fractions.end());
    ContinuationRequest request;
    request.start_problem = problem;
    request.target_problem = problem;
    for (std::size_t i = 0; i < 3; ++i) request.target_problem.feed.values[i] = 0.4 * x1[i] + 0.6 * x2[i];
    request.start_result = first;
    request.options = options(0.25, 0.25, 4, 2);
    epcsaft::core::Capture capture;
    const auto result = epcsaft::equilibrium::continue_equilibrium(model, request);
    record_bool("e4-e/path", "e4-e", "reactive two-liquid branch follows the feed path to the target", result.success, summary(result));
    const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
    for (std::size_t k = 0; k < result.states.size(); ++k) {
        const double q = result.path_values[k];
        const auto& state = result.states[k];
        double composition_error = 0.0, amount_error = 0.0, balance_error = 0.0, affinity = 0.0;
        for (std::size_t i = 0; i < 3; ++i)
            composition_error = std::max({composition_error, std::abs(state.mole_fractions[i] - x1[i]), std::abs(state.mole_fractions[3 + i] - x2[i])});
        amount_error = std::max(std::abs(state.amounts[0] - ((1.0 - q) * first.amounts[0] + 0.4 * q)),
                                std::abs(state.amounts[1] - ((1.0 - q) * first.amounts[1] + 0.6 * q)));
        std::vector<double> feed(3), totals(3, 0.0);
        for (std::size_t i = 0; i < 3; ++i) {
            feed[i] = (1.0 - q) * problem.feed.values[i] + q * request.target_problem.feed.values[i];
            totals[i] = state.amounts[0] * state.mole_fractions[i] + state.amounts[1] * state.mole_fractions[3 + i];
        }
        const auto lhs = epcsaft::equilibrium::matvec(compiled.B, totals), rhs = epcsaft::equilibrium::matvec(compiled.B, feed);
        for (std::size_t r = 0; r < lhs.size(); ++r) balance_error = std::max(balance_error, std::abs(lhs[r] - rhs[r]));
        std::array<double, 2> phase_affinity{};
        for (std::size_t row = 0; row < state.rows.size(); ++row)
            if (state.rows[row].kind == epcsaft::equilibrium::RowKind::Chemical)
                phase_affinity[state.rows[row].phase] += problem.reactions[0].stoichiometry[state.rows[row].index] * state.residuals[row];
        affinity = std::max(std::abs(phase_affinity[0]), std::abs(phase_affinity[1]));
        const double limit = std::max(state.convergence.criterion->value, state.convergence.floor.value_or(0.0));
        const bool verdict = composition_error <= 1.0e-8 && amount_error <= 1.0e-8 && balance_error <= 1.0e-12 && affinity <= limit;
        const auto id = std::string("e4-e/state/") + std::to_string(k);
        recording_support::compare(id.c_str(), "amount_error", "mol", amount_error, 0.0,
            recording_support::fields({recording_support::field("q", q), recording_support::field("amounts", state.amounts),
                recording_support::field("x", state.mole_fractions), recording_support::field("composition_error", composition_error),
                recording_support::field("balance_error", balance_error), recording_support::field("affinity", affinity),
                recording_support::field("limit", limit)}),
            "Relationship", "reactive-continuation", "reactive two-liquid amounts interpolate with unchanged compositions and zero affinity",
            "strict_absolute", 1.0e-8, 0.0, verdict, "e4-e", id.c_str(), "\\\\Delta N", "q", "\\\\Delta N\\\\leq10^{-8}", "q");
    }
    branch_rows("e4-e", "reactive", model, problem, result);
    attempt_rows("e4-e", capture.trace());
}

void check_e4_paths() {
    using namespace e4;
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    const auto opts = options(0.5, 0.5, 2, 1);
    const auto expect_invalid = [&](const char* id, ContinuationRequest request, const char* phrase) {
        request.options = opts;
        const auto result = epcsaft::equilibrium::continue_equilibrium(model, request);
        record_bool(id, "e4-paths", phrase, !result.success && result.stop == ContinuationStop::invalid_path
      && result.states.empty() && result.attempts == 0 && result.message.find(phrase) != std::string::npos, summary(result));
    };
    ContinuationRequest base;
    base.start_problem = equilibrium_support::bubble_at(180.0);
    base.target_problem = equilibrium_support::bubble_at(181.0);
    auto swapped = base; std::swap(swapped.target_problem.phases[0], swapped.target_problem.phases[1]);
    expect_invalid("e4-paths/phase-order", swapped, "identity differs");
    auto support = base; support.target_problem.phases[0].support = {0};
    expect_invalid("e4-paths/support", support, "support differs");
    auto reaction = base; reaction.target_problem.reactions = {{"A=B", {-1.0, 1.0}, 2.0, std::nullopt}};
    expect_invalid("e4-paths/reaction", reaction, "reaction count differs");
    auto law = base;
    law.start_problem.reactions = {{"A=B", {-1.0, 1.0}, 1.0, std::nullopt}};
    law.target_problem.reactions = {{"A=B", {-1.0, 1.0}, 2.0, std::nullopt}};
    expect_invalid("e4-paths/reaction-law", law, "reaction A=B differs");
    // Exact all-pinned endpoints pass pin admission but lose generic rank (E2-0 exact-sum record).
    ContinuationRequest pinned;
    pinned.start_problem = equilibrium_support::flash(2);
    pinned.start_problem.phases = {equilibrium_support::phase("L", PhaseKind::Liquid, AmountSpec::pinned(0.5)),
                                   equilibrium_support::phase("V", PhaseKind::Vapor, AmountSpec::pinned(0.5))};
    pinned.start_problem.temperature = FreeValue::variable(300.0);
    pinned.start_problem.pressure = FreeValue::variable(1.0e5);
    pinned.target_problem = pinned.start_problem;
    pinned.target_problem.feed.values = {1.2, 0.8};
    expect_invalid("e4-paths/rank-deficient-endpoint", pinned, "endpoint is not full rank");
    auto amount = base; amount.target_problem.phases[1].amount = AmountSpec::free();
    expect_invalid("e4-paths/amount-kind", amount, "amount declaration differs");
    auto flags = base; flags.target_problem.temperature = FreeValue::variable(181.0); flags.target_problem.pressure = FreeValue::fixed(1.0e5);
    expect_invalid("e4-paths/intensive-flags", flags, "temperature free/fixed declaration differs");
    auto tolerance = base; tolerance.target_problem.row_tolerance = 1.0e-8;
    expect_invalid("e4-paths/row-tolerance", tolerance, "row_tolerance differs");
    auto noop = base; noop.target_problem = noop.start_problem;
    expect_invalid("e4-paths/no-op", noop, "no fixed intensive or resolved feed differs");
    EquilibriumProblem liquid = equilibrium_support::flash(2);
    liquid.feed.values = {1.0, 0.0};
    ContinuationRequest zero_to_positive;
    zero_to_positive.start_problem = liquid;
    zero_to_positive.target_problem = liquid;
    zero_to_positive.target_problem.feed.values = {1.0, 1.0};
    expect_invalid("e4-paths/zero-to-positive", zero_to_positive, "admitted support of phase L differs");
    ContinuationRequest positive_to_zero;
    positive_to_zero.start_problem = liquid; positive_to_zero.start_problem.feed.values = {1.0, 1.0};
    positive_to_zero.target_problem = liquid;
    expect_invalid("e4-paths/positive-to-zero", positive_to_zero, "admitted support of phase L differs");
    // Real two-root admission falsifier: a selected liquid root cannot certify a vapor density.
    epcsaft::HelmholtzModel ethane(properties_support::ethane());
    const double T = 200.0, P = 217800.69494692178;
    const auto liquid_root = epcsaft::solve_density(ethane, T, P, {1.0}, epcsaft::DensityBranch::Liquid, std::nullopt);
    const auto vapor_root = epcsaft::solve_density(ethane, T, P, {1.0}, epcsaft::DensityBranch::Vapor, std::nullopt);
    const auto relabelled = epcsaft::solve_density(ethane, T, P, {1.0}, epcsaft::DensityBranch::Vapor, liquid_root.molar_density);
    const double displacement = std::abs(std::log(vapor_root.molar_density / liquid_root.molar_density));
    const double relabelled_displacement = std::abs(std::log(relabelled.molar_density / liquid_root.molar_density));
    const bool rejected = liquid_root.stable_root_count.value_or(0) == 2 && displacement > epcsaft::policy::DENSITY_ANCHOR_LOG_WINDOW.value
        && relabelled_displacement <= epcsaft::policy::DENSITY_ANCHOR_LOG_WINDOW.value
        && std::abs(std::log(vapor_root.molar_density / relabelled.molar_density)) > epcsaft::policy::DENSITY_ANCHOR_LOG_WINDOW.value;
    recording_support::compare_custom("e4-paths/two-root-admission", "log_displacement", "1", displacement,
        epcsaft::policy::DENSITY_ANCHOR_LOG_WINDOW.value,
        recording_support::fields({recording_support::field("T", T), recording_support::field("P", P),
            recording_support::field("rho_liquid", liquid_root.molar_density), recording_support::field("rho_vapor", vapor_root.molar_density),
            recording_support::field("rho_anchored_with_vapor_label", relabelled.molar_density),
            recording_support::field("relabelled_displacement", relabelled_displacement)}),
        "Relationship", "anchored-branch", "a corrected vapor density is outside the selected liquid root window even with a swapped label",
        "model > reference", rejected, "e4-paths", "two-root", "|\\\\ln(\\\\rho_V/\\\\rho_L)|", "T,P", "|\\\\ln(\\\\rho_V/\\\\rho_L)|>w", "case");
}

}  // namespace

int main(int argc, char** argv) {
    return recording_support::run(argc, argv, {
        {{"native.boundary.s3", "S3"}, check_s3},
        {{"native.boundary.s4", "S4"}, check_s4},
        {{"native.boundary.s5", "S5"}, check_s5},
        {{"native.boundary.s9", "S9"}, check_s9},
        {{"native.boundary.ionic-packing-derivative", "packing-derivative"}, check_ionic_packing_derivative},
        {{"native.boundary.s10", "S10"}, check_s10},
        {{"native.boundary.s14", "S14"}, check_s14},
        {{"native.boundary.s15c", "S15c"}, check_s15c},
        {{"native.boundary.e4-a", "E4-A"}, check_e4_a},
        {{"native.boundary.e4-b", "E4-B"}, check_e4_b},
        {{"native.boundary.e4-c", "E4-C"}, check_e4_c},
        {{"native.boundary.e4-d", "E4-D"}, check_e4_d},
        {{"native.boundary.e4-e", "E4-E"}, check_e4_e},
        {{"native.boundary.e4-paths", "E4-P"}, check_e4_paths},
    });
}
