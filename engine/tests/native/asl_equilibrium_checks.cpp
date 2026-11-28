#include "../../native/bindings/asl_adapter.cpp"
#include "equilibrium_support.hpp"

namespace {

std::string reactive_record() {
    // The same manufactured reactive LLE as e5-3, expressed in the ASL record's declared coordinates.
    const auto p = equilibrium_support::synthetic_reactive_ternary();
    epcsaft::ParameterDocument doc;
    doc.component_ids = {"A", "B", "C"};
    doc.base_choice = "pc-saft";
    doc.association_choice = "general-site";
    doc.electrolyte_choice = doc.permittivity_choice = "none";
    doc.site_component = p.site_component;
    doc.site_multiplicity = p.site_multiplicity;
    const auto entry = [&](const char* name, std::vector<std::size_t> indices, double value, const char* unit) {
        doc.entries.push_back({name, name, std::move(indices), value, unit});
    };
    for (std::size_t i = 0; i < 3; ++i) {
        entry("molar_mass", {i}, p.molar_mass_kg_per_mol[i], "kilogram / mole");
        entry("charge_number", {i}, 0.0, "elementary-charge");
        entry("segment_count", {i}, p.segment_count[i], "dimensionless");
        entry("sigma_constant", {i}, p.sigma_constant[i], "angstrom");
        entry("dispersion_energy_over_k", {i}, p.dispersion_energy_over_k[i], "kelvin");
        for (std::size_t j = i + 1; j < 3; ++j) entry("k_ij", {i, j}, p.k_ij[3 * i + j], "dimensionless");
    }
    entry("no_cross_association", {0, 1}, 0.0, "dimensionless");
    for (const std::size_t k : {0, 2}) {
        entry("association_energy_over_k", {k, k + 1}, p.association_energy_over_k[5 * k + 1], "kelvin");
        entry("association_volume", {k, k + 1}, p.association_volume[5 * k + 1], "dimensionless");
    }
    auto problem = equilibrium_support::synthetic_reactive_problem();
    problem.topology_declared = true;
    for (auto& phase : problem.phases) {
        phase.support = {0, 1, 2};
        phase.support_declared = true;
    }
    return eq::asl_record(doc, {}, {}, problem, {}, {{0, 0, 0}, problem.feed.values, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, 1, 0.1}, {});
}

struct Callback final {
    std::string record = reactive_record();
    std::array<double, 5> q{298.15, 1.0e5, 0.05, 0.475, 0.475};
    std::array<double, 5> gradient{};
    std::array<double, 15> hessian{};
    std::array<char, 5> fixed{1, 1, 0, 0, 0};

    double operator()(int order, int row = 0) {
        int map[]{0, 1, 2, 3, 4, -1, -2, -3};
        const char* ordinals[]{"0", "1", "2"};
        const char* strings[]{record.c_str(), ordinals[row / 3], ordinals[row % 3]};
        arglist args{};
        args.n = 8; args.nr = 5; args.at = map; args.ra = q.data(); args.sa = strings;
        args.dig = fixed.data();
        if (order >= 1) args.derivs = gradient.data();
        if (order >= 2) args.hes = hessian.data();
        const double result = equilibrium_observable(&args, 0);
        CHECK(args.Errmsg == nullptr && std::isfinite(result));
        return result;
    }
};

std::size_t spans(const epcsaft::core::Capture& trace, std::string_view name) {
    return std::count_if(trace.trace().spans.begin(), trace.trace().spans.end(), [&](const auto& span) { return span.name == name; });
}

std::vector<double> phase_component_amounts(const eq::EquilibriumResult& result) {
    auto amounts = result.mole_fractions;
    const std::size_t components = amounts.size() / result.phase_count;
    for (std::size_t i = 0; i < amounts.size(); ++i) amounts[i] *= result.amounts[i / components];
    return amounts;
}

void compare_amounts(const char* name, const EquilibriumContext& context, const eq::EquilibriumResult& warm,
                     const eq::EquilibriumResult& cold) {
    CHECK(warm.has_state() && cold.has_state());
    const auto left = phase_component_amounts(warm), right = phase_component_amounts(cold);
    const double tolerance = context.anchor.solver_options.value(epcsaft::policy::row)
        * epcsaft::policy::EQUILIBRIUM_INNER_TOLERANCE_FACTOR.value;
    CHECK(left.size() == right.size());
    for (std::size_t i = 0; i < left.size(); ++i) {
        recording_support::compare_custom(name + std::string("/") + std::to_string(i), "phase_component_amount", "mol",
            left[i], right[i], "inner solver tolerance at a 1 mol feed", std::abs(left[i] - right[i]) <= tolerance, "reactive-liquid-liquid-equilibrium");
        CHECK(std::abs(left[i] - right[i]) <= tolerance);
    }
}

eq::EquilibriumResult& selected(EquilibriumContext& context, const std::vector<double>& q) {
    return *select_state(context, q, "phase_component_amount").central;
}

void check_warm_start() {
    Callback callback;
    const auto context = prepare_equilibrium(callback.record.c_str(), "phase_component_amount");
    auto cold = prepare_equilibrium(callback.record.c_str(), "phase_component_amount");
    auto q = std::vector<double>(callback.q.begin(), callback.q.end());
    q[2] += 0.001;
    auto& first = selected(*context, q);
    CHECK(first.optimizer_state);
    q[2] += 0.001;
    const auto& warm = selected(*context, q);
    const auto& reference = selected(*cold, q);
    CHECK(warm.attempt == eq::Attempt::Warm);
    compare_amounts("reactive-lle/warm", *context, warm, reference);
}

void check_cold_fallback() {
    Callback callback;
    const auto context = prepare_equilibrium(callback.record.c_str(), "phase_component_amount");
    auto cold = prepare_equilibrium(callback.record.c_str(), "phase_component_amount");
    auto q = std::vector<double>(callback.q.begin(), callback.q.end());
    q[2] += 0.003;
    auto& first = selected(*context, q);
    CHECK(first.optimizer_state);
    auto one_liquid = context->anchor;
    one_liquid.phases.resize(1);
    const auto incompatible = eq::solve_equilibrium(*context->model, one_liquid);
    CHECK(incompatible.optimizer_state);
    first.optimizer_state = incompatible.optimizer_state;  // force the native dual-layout check to refuse
    q[2] += 0.001;
    const auto& reference = selected(*cold, q);
    epcsaft::core::Capture trace(true);
    const auto& result = selected(*context, q);
    CHECK(std::any_of(trace.trace().events.begin(), trace.trace().events.end(), [](const auto& event) {
        return std::string_view(event.name) == "warm_ineligible" && event.reason == "layout changed";
    }));
    CHECK(spans(trace, "solve_equilibrium") == 2);
    CHECK(result.attempt == eq::Attempt::Cold);
    compare_amounts("reactive-lle/cold-fallback", *context, result, reference);
}

void check_interleaved_cache() {
    Callback callback;
    const double value = callback(0);
    CHECK(callback(1) == value);
    CHECK(std::isfinite(callback.gradient[2]));
    const auto gradient = callback.gradient;
    CHECK(callback(2) == value);
    CHECK(callback.gradient == gradient);
    CHECK(std::isfinite(callback.hessian[5]));
    const double initial_feed = callback.q[2];
    for (int i = 1; i <= 4; ++i) {
        callback.q[2] = initial_feed + 0.001 * i;
        callback(0);
    }
    callback.q[2] = initial_feed;
    epcsaft::core::Capture trace(true);
    CHECK(callback(2) == value);
    CHECK(spans(trace, "solve_equilibrium") == 0);
    CHECK(spans(trace, "solved_state_actions") == 0);
}

}  // namespace

int main(int argc, char** argv) {
    return recording_support::run(argc, argv, {
        {{"native.regression.asl-reactive-lle-warm", "ASLReactiveLLEWarm"}, check_warm_start},
        {{"native.regression.asl-reactive-lle-cold-fallback", "ASLReactiveLLEColdFallback"}, check_cold_fallback},
        {{"native.regression.asl-reactive-lle-interleaved-cache", "ASLReactiveLLEInterleavedCache"}, check_interleaved_cache},
    });
}
