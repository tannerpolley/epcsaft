#define main original_properties_main
#include "properties_checks.cpp"
#undef main
#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium_support.hpp"
#include "site_input.hpp"
#include "fit_input.hpp"
#include <new>
#include <cstring>
inline thread_local long fail_after = -1;
inline thread_local unsigned injected = 0;
void* operator new(std::size_t size) {
    const bool eager_ipopt = epcsaft::core::ipopt_payload_call && !epcsaft::core::observation_depth;
    if (eager_ipopt) ++epcsaft::core::ipopt_eager_allocations;
    const bool endpoint = fail_after == -2 && epcsaft::core::observation_depth == 1
        && epcsaft::core::active_trace && !epcsaft::core::active_trace->diagnostics;
    const bool payload = fail_after == -2 - epcsaft::core::payload_kind;
    if ((eager_ipopt && fail_after == -5) || (payload && epcsaft::core::payload_kind) || endpoint
        || (epcsaft::core::observation_depth && fail_after >= 0 && fail_after-- == 0)) {
        fail_after = -1;
        ++injected;
        throw std::bad_alloc();
    }
    if (auto p = std::malloc(size ? size : 1)) return p;
    throw std::bad_alloc();
}
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
namespace e = epcsaft;
namespace c = epcsaft::core;
struct Outcome { std::vector<double> values; std::string cause; std::array<std::size_t, 10> work; };
Outcome execute(int which, long fault, bool capture_on = true) {
    const auto fit_request = which == 5 ? std::optional(recovery_request(manufactured(), 148.0, 0.0)) : std::nullopt;
    e::HelmholtzModel model(which == 4 ? properties_support::aqueous_electrolyte() : which == 2 ? site_input() : which == 0 ? properties_support::neutral_binary()
        : which == 1 ? p4_association_stress(2, 1e26) : properties_support::associating_binary());
    c::qualification_work = {}; c::payload_calls = {}; c::ipopt_eager_allocations = 0; injected = 0;
    std::optional<c::Capture> parent, capture;
    if (capture_on) { parent.emplace(true); capture.emplace(); }
    Outcome out;
    std::array<std::size_t, 3> driver{};
    fail_after = fault;
    try {
        if (which == 0) {
            const auto r = e::equilibrium::solve_equilibrium(model, equilibrium_support::flash(2));
            driver = {r.iterations, r.evaluations, r.hessian_calls};
            out.values = r.coordinates; out.values.insert(out.values.end(), r.residuals.begin(), r.residuals.end());
            out.values.insert(out.values.end(), r.classified_residuals.begin(), r.classified_residuals.end());
            out.cause = r.message;
        } else if (which == 5) {
            const auto r = e::regression::fit(*fit_request);
            CHECK(r.status == "converged" && r.usable && r.trial_failures == 0);
            driver = {r.iterations, r.residual_evaluations, r.jacobian_evaluations};
            out.values = r.physical; out.values.insert(out.values.end(), r.weighted_residuals.begin(), r.weighted_residuals.end());
            out.values.insert(out.values.end(), r.optimizer_jacobian.begin(), r.optimizer_jacobian.end());
            out.values.push_back(r.final_cost); out.cause = r.status + "/" + r.termination + "/" + r.covariance.status;
        } else if (which == 4) out.values = {model.evaluate(300, .001, {.98,.01,.01}, e::DerivativeOrder::Hessian).value};
        else if (which == 1) model.evaluate(300, .001, {1}, e::DerivativeOrder::Hessian);
        else if (which == 2) out.values = {model.state_at_pressure(45, 1e-4, {.5, .5}, e::DensityBranch::Liquid).molar_density};
        else {
            const auto r = e::reaction_references(model, 313.15, 1e5, {1, {1,0}, 1}, {{1,-1}}, {{}}, false, 1e30);
            out.values = {r[0].terminal_contraction.value(), r[0].solvent_density}; out.cause = r[0].message;
        }
    } catch (const std::exception& error) { out.cause = error.what(); }
    const auto w = model.counters();
    out.work = {(which == 0 || which == 5) ? driver[0] : w.forward, (which == 0 || which == 5) ? driver[1] : w.reverse, (which == 0 || which == 5) ? driver[2] : w.hessian_calls,
        which == 0 ? 0 : w.recordings, w.directional_paths, w.maximum_energy_order,
        c::qualification_work[0], c::qualification_work[1], c::qualification_work[2], c::qualification_work[3]};
    const auto t = capture ? capture->take() : c::Trace{};
    fail_after = -1;
    if (injected) CHECK(t.incomplete || parent->trace().incomplete || std::any_of(t.events.begin(), t.events.end(), [](const auto& v) { return v.reason == "std::bad_alloc"; }));
    if (capture_on && (fault == -4 || (fault == -3 && which == 4))) CHECK(injected == 1 && t.incomplete);
    if (capture_on && ((fault == -5 && which == 0) || (fault <= -6 && which == 5))) CHECK(injected == 1 && t.incomplete);
    if (!capture_on) CHECK(injected == 0 && !t.incomplete && c::payload_calls == (std::array<std::size_t, 6>{}));
    CHECK(c::ipopt_eager_allocations == 0);
    if (fault == -1) CHECK(!t.incomplete);
    for (const auto& v : t.events) CHECK(v.span < t.spans.size());
    if (which == 3 && fault == -1) for (const auto& v : t.events)
        if (std::strncmp(v.name, "anchor", 6) == 0) std::printf("anchor,%s,%.17g,%s\n", v.name, v.values[0], v.reason.c_str());
    if (fault == -2) for (const auto& v : t.events)
        if (std::string(v.name) == "solved_endpoint_jacobian") std::printf("endpoint-unavailable,%.17g,%s\n", v.values[2], v.reason.c_str());
    if (fault == -1 && !out.values.empty()) std::printf("returned,%d,%.17g\n", which, out.values.back());
    if (which == 5 && fault == -1) std::printf("fit-returned,%d,%.17g,%.17g,%.17g,%s,%zu,%zu,%zu\n",
        capture_on, out.values[0], out.values[1], out.values.back(), out.cause.c_str(), driver[0], driver[1], driver[2]);
    std::printf("allocation,%d,%ld,%u,%d,%zu,%zu", which, fault, injected, t.incomplete, t.spans.size(), t.events.size());
    for (auto value : out.work) std::printf(",%zu", value);
    std::printf(",%d", capture_on);
    for (const auto calls : c::payload_calls) std::printf(",%zu", calls);
    std::printf(",%zu\n", c::ipopt_eager_allocations);
    return out;
}
int main() {
    for (auto [sites, q] : std::array<std::pair<std::size_t,double>,2>{{{2,1e26},{128,1e27}}}) {
        e::HelmholtzModel model(p4_association_stress(sites, q), {}, p3_record(p3_forms()[0].second));
        c::Capture capture(true);
        const auto batch = model.density_actions(300, 1000, {1}, {e::PropertyObservable::ResidualEnthalpy}, {{.1,0,{0},{}}});
        const auto& events = capture.trace().events;
        auto event = std::find_if(events.begin(), events.end(), [](const auto& v) { return std::string(v.name) == "association_derivative_failure"; });
        CHECK(event != events.end() && event->values[1] == batch.diagnostics.site_rank && batch.diagnostics.site_factorizations == 1);
        CHECK(std::count_if(events.begin(), events.end(), [](const auto& v) { return std::string(v.name) == "association_derivative_failure"; }) == 1);
        if (std::isfinite(event->values[4])) CHECK(close(event->values[4], batch.diagnostics.site_pivot_ratio, 1e-12));
        std::printf("derivative,%zu,%zu,%.17g,%.17g,%.17g,%s\n", sites, batch.diagnostics.site_rank, event->values[2], event->values[3], event->values[4], event->reason.c_str());
        try { model.evaluate(300, .001, {1}, e::DerivativeOrder::Hessian); CHECK(false); }
        catch (const std::exception& error) { CHECK(std::string(error.what()).find("association Jacobian") != std::string::npos); }
    }
    for (auto [which, capture_on] : std::array<std::pair<int,bool>,9>{{{0,true},{1,true},{2,true},{3,true},{4,true},{4,false},{0,false},{5,true},{5,false}}}) {
        const auto baseline = execute(which, -1, capture_on);
        for (long fault : {-8,-7,-6,-5,-4,-3,-2,0,1,2,3,7,31,63,127,511,2047}) {
            if ((fault <= -6 && which != 5) || (fault == -5 && which != 0)) continue;
            if (!capture_on && fault >= -2) continue;
            const auto probe = execute(which, fault, capture_on);
            CHECK(probe.cause == baseline.cause && probe.work == baseline.work && probe.values.size() == baseline.values.size());
            CHECK(std::memcmp(probe.values.data(), baseline.values.data(), probe.values.size() * sizeof(double)) == 0);
        }
    }
    CHECK(c::active_trace == nullptr);
    return 0;
}
