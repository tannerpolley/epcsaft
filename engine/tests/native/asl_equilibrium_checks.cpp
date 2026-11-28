#include "../../native/bindings/asl_adapter.cpp"
#include "equilibrium_support.hpp"

namespace {

std::string reactive_record(bool warm_start = false) {
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
    return eq::asl_record(doc, {}, {}, problem, {}, {{0, 0, 0}, problem.feed.values, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
                                                    1, 0.1, warm_start}, {});
}

struct Callback final {
    std::string record = reactive_record();
    std::array<double, 5> q{298.15, 1.0e5, 0.05, 0.475, 0.475};
    std::array<double, 5> gradient{};
    std::array<double, 15> hessian{};
    std::array<char, 5> fixed{0, 1, 0, 0, 0};  // live T: at fixed T and P the reactive LLE has no degree of freedom, so the feed Hessian is zero

    double operator()(int order, int row = 0, bool enthalpy = false) {
        int map[]{0, 1, 2, 3, 4, -1, -2, -3};
        const char* ordinals[]{"0", "1", "2"};
        const char* strings[]{record.c_str(), ordinals[row / 3], ordinals[row % 3]};
        arglist args{};
        args.n = enthalpy ? 6 : 8; args.nr = 5; args.at = map; args.ra = q.data(); args.sa = strings;
        args.dig = fixed.data();
        if (order >= 1) args.derivs = gradient.data();
        if (order >= 2) args.hes = hessian.data();
        const double result = equilibrium_observable(&args, enthalpy ? 1 : 0);
        if (enthalpy) {
            CHECK(args.Errmsg != nullptr && std::isnan(result));
            CHECK(std::string_view(args.Errmsg).find("missing_thermochemistry:system.enthalpy") != std::string_view::npos);
        } else CHECK(args.Errmsg == nullptr && std::isfinite(result));
        return result;
    }
};

std::size_t spans(const epcsaft::core::Capture& trace, std::string_view name) {
    return std::count_if(trace.trace().spans.begin(), trace.trace().spans.end(), [&](const auto& span) { return span.name == name; });
}

void check_interleaved_cache() {
    Callback callback;
    callback.q[2] += 0.001;
    const double value = callback(0);
    CHECK(callback(1) == value);
    CHECK(std::isfinite(callback.gradient[2]));
    const auto gradient = callback.gradient;
    CHECK(callback(2) == value);
    CHECK(callback.gradient == gradient);
    CHECK(std::isfinite(callback.hessian[5]));
    const auto hessian = callback.hessian;
    const auto nonzero = [](const auto& values) { return std::any_of(values.begin(), values.end(), [](double v) { return v != 0.0; }); };
    CHECK(nonzero(gradient) && nonzero(hessian));
    const double initial_feed = callback.q[2];
    for (int i = 1; i <= 4; ++i) {
        callback.q[2] = initial_feed + 0.001 * i;
        callback(0);
    }
    callback.q[2] = initial_feed;
    callback.gradient.fill(0.0);  // the hit must rewrite every live derivative
    callback.hessian.fill(0.0);
    epcsaft::core::Capture trace(true);
    CHECK(callback(2) == value);
    CHECK(callback.gradient == gradient);
    CHECK(callback.hessian == hessian);
    CHECK(spans(trace, "solve_equilibrium") == 0);
    CHECK(spans(trace, "solved_state_actions") == 0);
}

void check_warm_start() {
    Callback cold, warm{reactive_record(true)};
    cold.q[2] += 0.006;
    warm.q = cold.q;
    warm(0);
    cold(0);
    cold.q[2] += 0.001;
    warm.q = cold.q;

    std::array<double, 6> cold_amounts{}, warm_amounts{};
    {
        epcsaft::core::Capture trace(true);
        cold(1);
        for (int i = 0; i < 6; ++i) cold_amounts[i] = cold(0, i);
        CHECK(spans(trace, "solve_equilibrium") == 1);
        CHECK(spans(trace, "equilibrium_warm_attempt") == 0);
    }
    {
        epcsaft::core::Capture trace(true);
        warm(1);
        for (int i = 0; i < 6; ++i) warm_amounts[i] = warm(0, i);
        CHECK(spans(trace, "solve_equilibrium") == 1);
        CHECK(spans(trace, "equilibrium_warm_attempt") == 1);
    }
    const auto scale = *std::max_element(cold_amounts.begin(), cold_amounts.end());
    for (std::size_t i = 0; i < cold_amounts.size(); ++i) {
        const bool agrees = std::abs(warm_amounts[i] - cold_amounts[i]) <= 1e-10 * scale;
        recording_support::compare_custom("reactive-lle/warm/" + std::to_string(i), "phase_component_amount", "mol",
            warm_amounts[i], cold_amounts[i], "1e-10 of cold amount scale", agrees, "reactive-liquid-liquid-equilibrium");
        CHECK(agrees);
    }
}

void check_cold_fallback() {
    const auto cold = prepare_equilibrium(reactive_record().c_str(), "phase_component_amount");
    const auto warm = prepare_equilibrium(reactive_record(true).c_str(), "phase_component_amount");
    std::vector<double> q{298.15, 1e5, 0.053, 0.475, 0.475};
    auto& seed = select_state(*warm, q, "phase_component_amount").central;
    CHECK(seed.optimizer_state);
    auto one_liquid = warm->anchor;
    one_liquid.phases.resize(1);
    const auto incompatible = eq::solve_equilibrium(*warm->model, one_liquid);
    CHECK(incompatible.optimizer_state);
    seed.optimizer_state = incompatible.optimizer_state;  // Force refusal by the native optimizer layout check.
    q[2] += 0.001;
    const auto& reference = select_state(*cold, q, "phase_component_amount").central;
    epcsaft::core::Capture trace(true);
    const auto& result = select_state(*warm, q, "phase_component_amount").central;
    CHECK(std::any_of(trace.trace().events.begin(), trace.trace().events.end(), [](const auto& event) {
        return std::string_view(event.name) == "warm_ineligible" && event.reason == "layout changed";
    }));
    CHECK(spans(trace, "solve_equilibrium") == 2);
    CHECK(result.attempt == eq::Attempt::Cold);
    CHECK(bitwise_equal(result.amounts, reference.amounts));
    CHECK(bitwise_equal(result.mole_fractions, reference.mole_fractions));
}

void check_refused_action_seed() {
    Callback warm{reactive_record(true)};
    warm.q[0] += 1.0;  // Different T isolates this sequence from other traces at the same record.
    warm.q[2] += 0.001;
    warm(1, 0, true);  // Equilibrium exists, but the requested enthalpy action lacks thermochemistry.
    {
        epcsaft::core::Capture trace(true);
        warm(0);
        CHECK(spans(trace, "solve_equilibrium") == 0);
    }
    warm.q[2] += 0.001;
    {
        epcsaft::core::Capture trace(true);
        warm(0);
        CHECK(spans(trace, "solve_equilibrium") == 1);
        CHECK(spans(trace, "equilibrium_warm_attempt") == 0);
    }
    warm.q[2] += 0.001;
    epcsaft::core::Capture trace(true);
    warm(0);
    CHECK(spans(trace, "equilibrium_warm_attempt") == 1);
}

}  // namespace

int main(int argc, char** argv) {
    return recording_support::run(argc, argv, {
        {{"native.regression.asl-reactive-lle-interleaved-cache", "ASLReactiveLLEInterleavedCache"}, check_interleaved_cache},
        {{"native.regression.asl-reactive-lle-warm-start", "ASLReactiveLLEWarmStart"}, check_warm_start},
        {{"native.regression.asl-reactive-lle-cold-fallback", "ASLReactiveLLEColdFallback"}, check_cold_fallback},
        {{"native.regression.asl-reactive-lle-refused-action-seed", "ASLReactiveLLERefusedActionSeed"}, check_refused_action_seed},
    });
}
