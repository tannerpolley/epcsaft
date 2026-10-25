#include "regression/observations.hpp"
#include "core/trace.hpp"
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <numeric>
namespace nb = nanobind;
using namespace nb::literals;
using namespace epcsaft;
using namespace epcsaft::equilibrium;
using namespace epcsaft::regression;
template <class... Args> nb::dict data(Args &&...args) {
    return nb::cast<nb::dict>(nb::module_::import_("builtins").attr("dict")(std::forward<Args>(args)...));
}
nb::list run(const HelmholtzModel &model, const std::vector<double> &formula, double T, double P) {
    nb::list records;
    std::vector<double> feed{formula.at(0), formula.at(1), formula.at(2), formula.at(3), formula.at(3)};
    for (int stage = 0; stage < 2; ++stage) {
        EquilibriumProblem problem{{}, FreeValue::fixed(T), FreeValue::fixed(P), {FeedKind::Amounts, feed}};
        problem.solver_options.set("EQUILIBRIUM_ROW_RESIDUAL", 1e-10);
        problem.phases = {{"aqueous", PhaseKind::Liquid, {0, 1, 2, 3, 4}, true, {}, {.80, .05, .05, .05, .05}}};
        if (stage)
            problem.phases.push_back(
                {"organic", PhaseKind::Liquid, {0, 1, 2, 3, 4}, true, {}, {.20, .10, .6998, .0001, .0001}});
        auto cell = compile_problem(model, problem);
        core::Capture capture(true);
        const auto s = solve_equilibrium(model, problem);
        nb::dict r = data("feed_mol"_a = feed, "feed_total"_a = cell.feed_total, "solved"_a = s,
                          "topology"_a = int(s.topology_event), "stop"_a = core::stop_name(s.convergence.stop));
        records.append(r);
        bool ok = s.success && s.has_state() && s.requested_tolerance_met == true &&
                  s.topology_event == TopologyEvent::None && maximum_abs(s.classified_residuals) <= 1e-10;
        r["admitted"] = ok;
        if (!ok) break;
        auto spans = [&]() {
            return std::count_if(capture.trace().spans.begin(), capture.trace().spans.end(),
                                 [](const auto &s) { return std::string(s.name) == "solve_equilibrium"; });
        };
        const auto before = spans();
        nb::list projections, phases;
        r["projections"] = projections;
        r["phases"] = phases;
        auto project = [&](std::size_t phase, std::vector<std::size_t> species, std::vector<std::size_t> denominator) {
            Observation o;
            o.problem = cell.problem;
            o.phase = o.denominator_phase = phase;
            o.species = species;
            o.denominator_species = denominator;
            o.kind = denominator.empty() ? ObservationKind::PhaseMoleFraction : ObservationKind::PhaseAmountRatio;
            Prediction p;
            p.compiled = cell;
            p.solved = s;
            p.convergence = s.convergence;
            p.requested_tolerance_met = s.requested_tolerance_met;
            p.message = s.message;
            p.solver_options = cell.problem.solver_options.resolved(model.solver_options());
            validate(o, model, &cell);
            project_prediction(model, o, p);
            const bool identical = p.solved->coordinates == s.coordinates && p.solved->amounts == s.amounts &&
                                   p.solved->mole_fractions == s.mole_fractions &&
                                   p.solved->molar_densities == s.molar_densities;
            projections.append(data("phase"_a = phase, "species"_a = species, "denominator"_a = denominator,
                                    "status"_a = observation_status_name(p.status), "value"_a = p.value.value_or(NAN),
                                    "identical"_a = identical, "message"_a = p.message));
            ok = ok && identical && p.status == ObservationStatus::Available && p.value && std::isfinite(*p.value);
            return p.value.value_or(NAN);
        };
        std::vector<double> balance = feed;
        double distance = 0;
        for (std::size_t phase = 0; phase < s.phase_count; ++phase) {
            std::vector<double> x(s.mole_fractions.begin() + phase * 5, s.mole_fractions.begin() + (phase + 1) * 5);
            const auto state = model.state(T, s.molar_densities[phase], x);
            const double amount = s.amounts[phase], pressure = std::abs(state.pressure / P - 1),
                         sum = std::accumulate(x.begin(), x.end(), 0.0);
            ok = ok && amount > 0 && std::isfinite(amount) && state.pressure_density_derivative > 0 &&
                 std::isfinite(state.pressure_density_derivative) && pressure <= 1e-8 && std::abs(sum - 1) <= 1e-10 &&
                 std::abs(amount * (x[3] - x[4])) <= 1e-10 && *std::min_element(x.begin(), x.end()) >= 0;
            nb::dict ph = data("pressure_relative"_a = pressure, "slope"_a = state.pressure_density_derivative,
                               "normalization"_a = sum - 1, "charge_mol"_a = amount * (x[3] - x[4]));
            phases.append(ph);
            for (std::size_t j = 0; j < 5; ++j) {
                balance[j] -= amount * x[j];
                if (stage) distance = std::max(distance, std::abs(s.mole_fractions[j] - s.mole_fractions[5 + j]));
            }
            for (int kind = 0; kind < 2; ++kind) {
                std::vector<std::size_t> denominator{0, 1, 2};
                if (!kind) denominator.push_back(3);
                const double fraction = project(phase, denominator, {});
                ok = ok && fraction > 0 && std::isfinite(fraction);
                std::vector<double> values;
                for (std::size_t j = 0; j < (kind ? 3 : 5); ++j) values.push_back(project(phase, {j}, denominator));
                ph[kind ? "salt_free" : "formula"] = values;
                if (!kind) {
                    ok = ok && std::abs(std::accumulate(values.begin(), values.begin() + 4, 0.0) - 1) <= 1e-10 &&
                         std::abs(values[3] - values[4]) <= 1e-10;
                    if (!stage)
                        for (std::size_t j = 0; j < 4; ++j) ok = ok && std::abs(values[j] - formula[j]) <= 1e-10;
                }
            }
        }
        ok = ok && maximum_abs(balance) <= 1e-10 && std::abs(feed[3] - feed[4]) <= 1e-10;
        if (stage)
            ok = ok && std::abs(cell.feed_total - 1) <= 1e-10 && distance >= 1e-3 &&
                 *std::min_element(s.amounts.begin(), s.amounts.end()) / cell.feed_total >= 1e-4;
        ok = ok && before > 0 && spans() == before && !capture.trace().incomplete;
        r.attr("update")(data("admitted"_a = ok, "balance_defects_mol"_a = balance, "phase_distance"_a = distance,
                              "solve_spans_before"_a = before, "solve_spans_after"_a = spans(),
                              "trace_incomplete"_a = capture.trace().incomplete));
        if (!ok) break;
        feed = s.mole_fractions;
    }
    return records;
}
NB_MODULE(figure2, m) {
    m.def("run", &run);
}
