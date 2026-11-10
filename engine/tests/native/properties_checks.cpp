#include "recording_support.hpp"
#include "equilibrium_support.hpp"
#include <condition_variable>
#include <future>
#include <thread>

// The case registry. Each case function is defined in one properties_*_checks.cpp unit.
namespace properties_checks {
void concurrent_first_discovery();
void residual_cv_temperature_extraction();
void born_parameter_columns();
void reference_f8();
void reference_f9();
void reference_limits();
void residual_properties_and_caloric_ladders();
void residual_internal_energy_identity_ladder();
void derivatives_and_limits();
void finite_difference_gradient();
void hessian_reference();
void hessian_fd_and_pressure_slope();
void same_sign_pair_is_inert();
void parameter_derivative_and_helpers();
void core_c2_parameter_coordinates();
void density_and_reference();
void properties_p2();
void properties_thermochemistry_p3();
void properties_p4();
void cache_invariance();
void modular_permittivity_actions();
void electrolyte_parameter_actions();
void energy_e5_ladder();
void energy_action_ladders();
void neutral_property_pressure_hessian();
void log_fugacity_actions();
void liquid_reference_branch_evidence();
void joule_thomson_identities();
void temperature_laws();
void cation_mixing();
void pathway_weight();
}  // namespace properties_checks

int main(int argc, char** argv) {
    using namespace properties_checks;
    return recording_support::run(argc, argv, {
        {{"native.properties.concurrent-first-discovery", "ConcurrentFirstDiscovery"}, concurrent_first_discovery},
        {{"native.properties.residual-cv-temperature-extraction", "ResidualCvTemperature"}, residual_cv_temperature_extraction},
        {{"native.properties.born-parameter", "BornParameter"}, born_parameter_columns},
        {{"native.properties.reference-f8", "F8"}, reference_f8, "F8"},
        {{"native.properties.reference-f9", "F9"}, reference_f9, "F9"},
        {{"native.properties.reference-limits", "ReferenceLimits"}, reference_limits, "ReferenceLimits"},
        {{"native.properties.p1", "P1"}, residual_properties_and_caloric_ladders, "P1"},
        {{"native.properties.f1", "F1"}, residual_internal_energy_identity_ladder, "F1"},
        {{"native.properties.m2", "M2"}, derivatives_and_limits},
        {{"native.properties.m3", "M3"}, finite_difference_gradient},
        {{"native.properties.m4", "M4"}, hessian_reference},
        {{"native.properties.m5", "M5"}, hessian_fd_and_pressure_slope},
        {{"native.properties.m6", "M6"}, same_sign_pair_is_inert},
        {{"native.properties.m7", "M7"}, parameter_derivative_and_helpers},
        {{"native.properties.c2", "C2"}, core_c2_parameter_coordinates},
        {{"native.properties.m8", "M8"}, density_and_reference},
        {{"native.properties.p2", "P2"}, properties_p2, "P2"},
        {{"native.properties.p3", "P3"}, properties_thermochemistry_p3, "P3"},
        {{"native.properties.p4", "P4"}, properties_p4, "P4"},
        {{"native.properties.m9", "M9"}, cache_invariance},
        {{"native.properties.issue26", "Issue26"}, modular_permittivity_actions},
        {{"native.properties.electrolyte-coordinates", "ElectrolyteCoordinates"}, electrolyte_parameter_actions},
        {{"native.properties.e5-ladder", "C4"}, energy_e5_ladder},
        {{"native.properties.issue26e", "Issue26E"}, energy_action_ladders},
        {{"native.properties.issue10", "Issue10"}, neutral_property_pressure_hessian},
        {{"native.properties.log-fugacity", "LogFugacity"}, log_fugacity_actions},
        {{"native.properties.cation-mixing", "CationMixing"}, cation_mixing},
        {{"native.properties.temperature-laws", "TemperatureLaws"}, temperature_laws},
        {{"native.properties.pathway-weight", "PathwayWeight"}, pathway_weight},
        {{"native.properties.issue423", "Issue423"}, liquid_reference_branch_evidence},
        {{"native.properties.joule-thomson", "JouleThomson"}, joule_thomson_identities},
    });
}

namespace properties_checks {
void concurrent_first_discovery() {
    using namespace epcsaft;
    using namespace epcsaft::equilibrium;
    const auto problem = [](bool ternary) {
        auto p = equilibrium_support::flash(ternary ? 3 : 2);
        p.temperature = FreeValue::fixed(190);
        p.pressure = FreeValue::fixed(1e6);
        p.feed.values = ternary ? std::vector<double>{.3, .3, .4} : std::vector<double>{.5, .5};
        p.phases.clear();
        p.topology_declared = false;
        p.discovery = StabilityOptions{};
        return p;
    };
    const auto parameters = [](bool ternary) {
        return ternary ? properties_support::neutral_ternary() : properties_support::neutral_binary();
    };
    using Work = decltype(EquilibriumResult{}.discovery->work);
    std::array<Work, 4> observed;
    std::mutex mutex;
    std::condition_variable ready;
    int arrived = 0;
    std::vector<std::future<void>> tasks;
    for (std::size_t worker = 0; worker < observed.size(); ++worker) {
        tasks.push_back(std::async(std::launch::async, [&, worker] {
            HelmholtzModel model(parameters(worker % 2));
            const auto p = problem(worker % 2);
            {
                std::unique_lock<std::mutex> lock(mutex);
                ++arrived;
                ready.notify_all();
                ready.wait(lock, [&] { return arrived == 4; });
            }
            const auto first = discover_phases(model, p);
            if (!first.success) throw std::runtime_error("cold discovery failed");
            observed[worker] = first.discovery->work;
            for (int repeat = 0; repeat < 20; ++repeat) {
                const auto result = discover_phases(model, p);
                if (!result.success || result.discovery->work != observed[worker])
                    throw std::runtime_error("concurrent discovery work differs");
            }
        }));
    }
    for (auto& task : tasks) task.get();
    for (std::size_t worker = 0; worker < observed.size(); ++worker) {
        HelmholtzModel model(parameters(worker % 2));
        const auto serial = discover_phases(model, problem(worker % 2));
        if (!serial.success || serial.discovery->work != observed[worker])
            throw std::runtime_error("serial discovery work differs");
        for (const auto& [key, value] : serial.discovery->work)
            std::cout << worker << ',' << key << ',' << value << '\n';
    }

    const auto evaluate = [] {
        HelmholtzModel model(properties_support::neutral_binary());
        const auto value = model.evaluate(190., .001, {.5, .5}, DerivativeOrder::Hessian);
        if (!std::isfinite(value.value) || !std::all_of(value.hessian.begin(), value.hessian.end(),
                [](double entry) { return std::isfinite(entry); }))
            throw std::runtime_error("nonfinite thread-exit Hessian");
        return value.value;
    };
    const double reference = evaluate();
    struct Cleanup final {
        std::mutex& mutex;
        std::condition_variable& ready;
        bool& exiting;
        bool& other_ready;
        bool& go;
        const decltype(evaluate)& run;
        double reference;
        ~Cleanup() {
            {
                std::unique_lock<std::mutex> lock(mutex);
                exiting = true;
                ready.notify_all();
                ready.wait(lock, [&] { return other_ready; });
                go = true;
                ready.notify_all();
            }
            for (int repeat = 0; repeat < 100; ++repeat)
                if (run() != reference) std::abort();
        }
    };
    {
        bool exiting = false, other_ready = false, go = false;
        std::thread first([&] {
            thread_local Cleanup cleanup{mutex, ready, exiting, other_ready, go, evaluate, reference};
            (void)cleanup;
            if (evaluate() != reference) std::abort();
        });
        {
            std::unique_lock<std::mutex> lock(mutex);
            ready.wait(lock, [&] { return exiting; });
        }
        std::thread second([&] {
            if (evaluate() != reference) std::abort();
            {
                std::unique_lock<std::mutex> lock(mutex);
                other_ready = true;
                ready.notify_all();
                ready.wait(lock, [&] { return go; });
            }
            for (int repeat = 0; repeat < 100; ++repeat)
                if (evaluate() != reference) std::abort();
        });
        first.join();
        second.join();
    }
    // More than 48 successive callers also falsify process-lifetime slot retention.
    for (int caller = 0; caller < 50; ++caller) {
        std::thread worker([&] { if (evaluate() != reference) std::abort(); });
        worker.join();
    }
    std::cout << "thread-exit: overlapping teardown and 50 recycled callers matched serial values\n";
}
}  // namespace properties_checks
