// Disposable native code verification: include the unchanged final owner to reach Evaluator.
#include "../../../engine/native/regression/fit.cpp"
#include "../../../engine/tests/native/equilibrium_support.hpp"
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

using namespace epcsaft;
using namespace epcsaft::regression;

int main(int argc, char** argv) {
    assert(argc == 2);
    std::ofstream rows(std::string(argv[1]) + "/residual-promotion-values.csv");
    rows << std::setprecision(17)
         << "stage,epsilon_K,row,T_K,prediction_Pa,residual,weighted_residual,cost,scalar_reference_Pa,prediction_difference_Pa,physical_derivative_Pa_K,optimizer_jacobian\n";
    FitRequest request;
    request.parameters = properties_support::neutral_binary();
    request.parameters.dispersion_energy_over_k[0] = 148.0;
    request.coordinates.coordinates = {{{{ParameterFamily::DispersionEnergyOverK, {0}}}, 148.0, 1.0e308, 140.0, 170.0}};
    // Difference in units of 1 Pa is an admitted public residual transform. Only the coordinate
    // scale is extreme: exact pressure derivatives are finite but their optimizer chain overflows.
    for (const double T : {170.0, 180.0}) {
        Observation row;
        row.kind = ObservationKind::EquilibriumPressure;
        row.problem = equilibrium_support::boundary();
        row.problem.temperature = equilibrium::FreeValue::fixed(T);
        row.observed = 1.0e6;
        row.form = ResidualForm::Difference;
        row.scale = 1.0;
        request.observations.push_back(row);
    }
    request.weights = {1.0, 0.25};
    validate(request.coordinates, request.parameters);
    const HelmholtzModel admitted(request.parameters, request.coordinates.active_parameters());
    for (const auto& row : request.observations) validate(row, admitted);
    const auto groups = state_groups(request.observations);
    const std::vector<double> sqrt_weights{1.0, 0.5};
    FitResult result;
    ObservationEvaluator::Impl evaluator(request, result, groups, sqrt_weights);
    const auto point = request.coordinates.start(request.parameters);
    assert(point == std::vector<double>{0.0});
    core::Capture capture;
    const auto count = [&](std::size_t begin, const char* name) {
        std::size_t total = 0;
        for (std::size_t i = begin; i < capture.trace().spans.size(); ++i)
            if (std::string(capture.trace().spans[i].name) == name) ++total;
        return total;
    };
    std::vector<std::size_t> solves, compilations;
    const auto evaluate = [&](const std::vector<double>& u, bool jacobian, ObservationEvaluation& out) {
        const auto begin = capture.trace().spans.size();
        const bool ok = evaluator.evaluate(u, jacobian, out);
        solves.push_back(count(begin, "solve_equilibrium"));
        compilations.push_back(count(begin, "compile_problem"));
        return ok;
    };
    ObservationEvaluation before, failed, after, changed;
    assert(evaluate(point, false, before));
    assert(finite(before.predictions) && finite(before.residuals) && finite(before.weighted) && std::isfinite(before.cost));
    assert(!evaluate(point, true, failed));
    assert(result.failures.size() == 1);
    const auto& failure = result.failures.front();
    assert(failure.observation == 0 && failure.status == ObservationStatus::Nonfinite && failure.jacobian
           && failure.message == "optimizer Jacobian not finite");
    assert(finite(failed.gradient) && std::isfinite(failed.cost) && !std::isfinite(failed.optimizer_gradient[0]));
    assert(evaluate(point, false, after));
    assert(after.predictions == before.predictions && after.residuals == before.residuals
           && after.weighted == before.weighted && after.cost == before.cost);
    std::vector<double> changed_point{1.0e-308};
    assert(request.coordinates.within_bounds(changed_point));
    assert(evaluate(changed_point, false, changed));
    const auto record = [&](const char* stage, const std::vector<double>& u, const ObservationEvaluation& values) {
        const auto physical = request.coordinates.physical(u);
        const HelmholtzModel scalar(request.coordinates.apply(request.parameters, u));
        for (std::size_t i = 0; i < request.observations.size(); ++i) {
            const auto independent = predict(scalar, request.observations[i]);
            assert(independent.status == ObservationStatus::Available && independent.value && independent.residual);
            const double difference = values.predictions[i] - *independent.value;
            assert(std::abs(difference) <= 1.0e-10 * std::abs(*independent.value));
            rows << stage << ',' << physical[0] << ',' << i << ',' << request.observations[i].problem.temperature.value
                 << ',' << values.predictions[i] << ',' << values.residuals[i] << ',' << values.weighted[i]
                 << ',' << values.cost << ',' << *independent.value << ',' << difference << ',';
            if (!values.gradient.empty()) rows << values.gradient[i];
            rows << ',';
            if (!values.optimizer_gradient.empty() && (std::string(stage) != "jacobian_attempt" || i <= failure.observation))
                rows << values.optimizer_gradient[i];
            rows << '\n';
        }
    };
    record("residual_before", point, before);
    record("jacobian_attempt", point, failed);
    record("residual_after", point, after);
    record("changed_point", changed_point, changed);
    std::ofstream diagnostic(std::string(argv[1]) + "/residual-promotion-diagnostic.json");
    diagnostic << std::setprecision(17) << "{\"verification\":\"native code\",\"scale\":1e308,\"origin_K\":148,\"physical_bounds_K\":[140,170],"
               << "\"failure_row\":" << failure.observation << ",\"failure\":\"" << observation_status_name(failure.status)
               << "\",\"message\":\"" << failure.message << "\",\"jacobian_requested\":true,\"prediction_residual_cost_identical\":true,"
               << "\"central_solve_counts\":[" << solves[0] << ',' << solves[1] << ',' << solves[2] << ',' << solves[3]
               << "],\"compile_counts\":[" << compilations[0] << ',' << compilations[1] << ',' << compilations[2] << ',' << compilations[3]
               << "],\"count_stage_order\":[\"residual_before\",\"jacobian_attempt\",\"residual_after\",\"changed_point\"]}\n";
    std::cout << std::setprecision(17) << "before_cost=" << before.cost << " after_cost=" << after.cost
              << " changed_cost=" << changed.cost << " failure_row=" << failure.observation << " failure="
              << observation_status_name(failure.status) << " message=" << failure.message << '\n';
}
