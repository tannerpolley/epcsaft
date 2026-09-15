#include "properties/helmholtz_model.hpp"
#include "properties_support.hpp"
#include "recording_support.hpp"

#include "core/trace.hpp"
#include "contributions/association.hpp"
#include "contributions/helmholtz.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

namespace {

double median(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}

template <typename Work>
double timed(Work&& work) {
    const auto start = std::chrono::steady_clock::now();
    work();
    const auto stop = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::micro>(stop - start).count();
}

void print_evaluation_case(
    const char* label,
    epcsaft::HelmholtzModel& model,
    double temperature,
    const std::vector<double>& amounts,
    double volume,
    epcsaft::DerivativeOrder order
) {
    for (int warmup = 0; warmup < 3; ++warmup) {
        model.evaluate(temperature, volume, amounts, order);
    }
    std::vector<double> times;
    for (int sample = 0; sample < 20; ++sample) {
        model.reset_counters();
        times.push_back(timed([&] {
            model.evaluate(temperature, volume, amounts, order);
        }));
    }
    const auto counts = model.counters();
    std::printf("%s wall_us=%.3f sweeps=%zu hessian_calls=%zu recordings=%zu\n",
                label, median(times), counts.forward + counts.reverse,
                counts.hessian_calls, counts.recordings);
}

void print_density_case(
    const char* label,
    epcsaft::HelmholtzModel& model,
    double temperature,
    double pressure,
    const std::vector<double>& composition,
    bool anchored,
    epcsaft::DensityBranch branch
) {
    try {
        for (int warmup = 0; warmup < 3; ++warmup) {
            const auto root = epcsaft::solve_density(
                model, temperature, pressure, composition,
                branch, std::nullopt
            );
            if (anchored) {
                epcsaft::solve_density(
                    model, temperature, pressure, composition,
                    branch, root.molar_density
                );
            }
        }
    } catch (const std::exception&) {
        std::printf("%s cold=not measured\n", label);
        return;
    }
    std::vector<double> times;
    epcsaft::DensityRoot last;
    for (int sample = 0; sample < 20; ++sample) {
        model.reset_counters();
        times.push_back(timed([&] {
            last = epcsaft::solve_density(
                model, temperature, pressure, composition,
                branch, std::nullopt
            );
            if (anchored) {
                last = epcsaft::solve_density(
                    model, temperature, pressure, composition,
                    branch, last.molar_density
                );
            }
        }));
    }
    const auto counts = model.counters();
    std::printf("%s wall_us=%.3f sweeps=%zu hessian_calls=%zu recordings=%zu "
                "density_iterations=%zu\n", label, median(times),
                counts.forward + counts.reverse, counts.hessian_calls,
                counts.recordings, last.iterations);
    std::printf("%s convergence=%s\n", label,
                epcsaft::core::describe(last.convergence).c_str());
}

void print_reference_case(
    epcsaft::HelmholtzModel& model,
    double temperature,
    double pressure
) {
    const std::vector<double> solvent = {1.0, 0.0, 0.0};
    const std::array<std::size_t, 2> ions = {1, 2};
    const std::array<int, 2> stoichiometry = {1, 1};
    try {
        for (int warmup = 0; warmup < 3; ++warmup) {
            epcsaft::infinite_dilution_reference(
                model, temperature, pressure, solvent, ions, stoichiometry
            );
        }
    } catch (const std::exception&) {
        std::printf("F infinite_dilution_reference=not measured\n");
        return;
    }
    std::vector<double> times;
    for (int sample = 0; sample < 20; ++sample) {
        model.reset_counters();
        times.push_back(timed([&] {
            epcsaft::infinite_dilution_reference(
                model, temperature, pressure, solvent, ions, stoichiometry
            );
        }));
    }
    const auto counts = model.counters();
    std::printf("F wall_us=%.3f sweeps=%zu hessian_calls=%zu recordings=%zu\n",
                median(times), counts.forward + counts.reverse,
                counts.hessian_calls, counts.recordings);
}

void print_trace_case(epcsaft::HelmholtzModel& model) {
    const std::vector<double> composition = {
        0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0
    };
    double uncaptured_wall = 0.0;
    double captured_wall = 0.0;
    const auto solve = [&] {
        return model.state_at_pressure(
            313.15, 1.0e5, composition, epcsaft::DensityBranch::Liquid
        );
    };
    try {
        uncaptured_wall = timed(solve);
    } catch (const std::exception& error) {
        std::printf("G uncaptured_error=%s\n", error.what());
    }
    epcsaft::core::Capture capture;
    try {
        captured_wall = timed(solve);
    } catch (const std::exception& error) {
        std::printf("G captured_error=%s\n", error.what());
    }
    const auto trace = capture.take();
    std::size_t values_bytes = 0;
    for (const auto& event : trace.events) {
        values_bytes += event.values.size() * sizeof(double);
    }
    std::printf("G uncaptured_wall_us=%.3f captured_wall_us=%.3f spans=%zu events=%zu values_bytes=%zu\n",
                uncaptured_wall, captured_wall, trace.spans.size(),
                trace.events.size(), values_bytes);
}

}  // namespace

void run_floor() {
    std::printf("compiler=%s\n", __VERSION__);
    std::ifstream cpuinfo("/proc/cpuinfo");
    std::string line;
    while (std::getline(cpuinfo, line)) {
        if (line.rfind("model name", 0) == 0) {
            std::printf("cpu=%s\n", line.c_str());
            break;
        }
    }
    const char* threads = std::getenv("OMP_NUM_THREADS");
    std::printf("OMP_NUM_THREADS=%s\n", threads == nullptr ? "unset" : threads);
    epcsaft::HelmholtzModel model(
        properties_support::nine_species(),
        {{epcsaft::ParameterFamily::DispersionEnergyOverK, {0}}}
    );
    print_evaluation_case("A", model, 313.15,
               {0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0},
               1.0 / 3.0e4, epcsaft::DerivativeOrder::Hessian);
    const auto assoc = properties_support::associating_binary();
    epcsaft::HelmholtzModel associating(assoc);
    print_evaluation_case("C", associating, 376.88, {0.1385, 0.8615},
                          1.0 / 1.0e4, epcsaft::DerivativeOrder::Hessian);
    auto association_parameters = associating.parameters();
    const auto sigma = epcsaft::segment_sigma(
        376.88, association_parameters
    );
    const auto diameter = epcsaft::segment_diameter(
        376.88, association_parameters, sigma
    );
    epcsaft::core::Convergence association_record;
    epcsaft::solve_site_fractions(
        376.88, 1.0e4, {0.1385, 0.8615}, association_parameters,
        sigma, diameter, &association_record
    );
    std::printf("C association=%s\n",
                epcsaft::core::describe(association_record).c_str());
    epcsaft::HelmholtzModel density_model(properties_support::neutral_binary());
    print_density_case("E", density_model, 300.0, 1.0e5, {0.4, 0.6}, false,
                       epcsaft::DensityBranch::Vapor);
    std::printf("B reference tensor multiplier=not measured\n");
    epcsaft::HelmholtzModel nine_density(properties_support::nine_species());
    print_density_case("D", nine_density, 313.15, 1.0e5,
                       {0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0}, true,
                       epcsaft::DensityBranch::Liquid);
    print_trace_case(nine_density);
    epcsaft::HelmholtzModel aqueous(properties_support::aqueous_electrolyte());
    print_reference_case(aqueous, 298.15, 1.0e5);
}

int main(int argc, char** argv) {
    return recording_support::run(argc, argv, {
        {{"native.properties_floor", "Model_floor_measurements", "Tooling", "Properties", "Explicit", "Outcome_only", nullptr, "P5", nullptr, "engine/tests/native/properties_floor.cpp", "properties_floor", "run"}, run_floor},
    });
}
