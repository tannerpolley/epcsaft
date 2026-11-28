#include "properties/helmholtz_model.hpp"
#include "surface/parameter_builder.hpp"
#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium_record.hpp"

#include <charconv>
#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <list>
#include <memory>
#include <sstream>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>

// Use the process stderr for retained diagnostics, including direct Pyomo calls without AE stdio hooks.
#define NO_STDIO1
#define No_AE_redefs
#include "../../build/subprojects/packagefiles/asl/funcadd.h"

namespace {

using epcsaft::HelmholtzModel;
using epcsaft::equilibrium::model_at;

thread_local std::string callback_error;

constexpr double invalid_value = std::numeric_limits<double>::quiet_NaN();

// ASL reads no partial with respect to a real argument it marks in `dig` (a constant of the problem); those stay untouched.
bool live(const arglist* arguments, int i) { return arguments->dig == nullptr || arguments->dig[i] == 0; }

void invalidate(arglist* arguments) {
    if (arguments == nullptr) return;
    for (int j = 0, packed = 0; j < arguments->nr; ++j) {
        if (arguments->derivs != nullptr && live(arguments, j)) arguments->derivs[j] = invalid_value;
        for (int i = 0; i <= j; ++i, ++packed)
            if (arguments->hes != nullptr && live(arguments, i) && live(arguments, j)) arguments->hes[packed] = invalid_value;
    }
}

[[noreturn]] void reject(const std::string& message) {
    throw std::invalid_argument(message);
}

const char* string_argument(arglist* arguments, int position) {
    if (arguments == nullptr) reject("missing arguments");
    const int index = arguments->at[position];
    if (index >= 0) reject("invalid string argument map");
    const char* value = arguments->sa[-index - 1];
    if (value == nullptr) reject("missing string argument");
    return value;
}

// Real arguments first, then `strings` string arguments, in declaration order; returns the real count.
int argument_count(arglist* arguments, int strings, const char* message) {
    constexpr int maximum = 2 + 2 * epcsaft::equilibrium::asl_record_vector_limit;
    if (arguments == nullptr || arguments->nr < 1 || arguments->nr > maximum || arguments->n != arguments->nr + strings
        || arguments->at == nullptr || arguments->ra == nullptr || arguments->sa == nullptr)
        reject(message);
    const int n = arguments->nr;
    for (int i = 0; i < arguments->n; ++i)
        if (arguments->at[i] != (i < n ? i : n - 1 - i)) reject("invalid argument map");
    return n;
}

std::size_t selector(arglist* arguments, int position, std::size_t maximum) {
    const std::string s(string_argument(arguments, position));
    std::size_t v = 0;
    const auto parsed = std::from_chars(s.data(), s.data() + s.size(), v);
    if (s.empty() || parsed.ec != std::errc{} || parsed.ptr != s.data() + s.size() || v >= maximum)
        reject("invalid observable selector");
    return v;
}

// One context per distinct record content, prepared once and retained for the process lifetime.
template <class Context> struct Prepared final {
    std::mutex mutex;
    std::vector<std::unique_ptr<Context>> contexts;
    template <class Prepare> Context& find(std::string_view record, const char* kind, Prepare prepare) {
        for (const auto& item : contexts) if (item->record == record) return *item;
        auto& result = *contexts.emplace_back(prepare());
        std::fprintf(stderr, "epcsaft_%s prepared record %zu\n", kind, contexts.size());
        return result;
    }
};

namespace eq = epcsaft::equilibrium;
using epcsaft::core::bitwise_equal;
using ActionStatus = eq::SolvedStateActionStatus;

class EquilibriumFailure final : public std::runtime_error {
public:
    EquilibriumFailure(ActionStatus status, const std::string& observable, const std::string& detail = {})
        : std::runtime_error(std::string("epcsaft_equilibrium:")
              + eq::solved_state_status_name(status) + ":" + observable + (detail.empty() ? "" : ": " + detail)) {}
};

// Record admission reconstructs native inputs; the Surface builder owns parameter rules.
class RecordReader final {
public:
    explicit RecordReader(const char* raw) : text_(raw, strnlen(raw, max_record + 1)) {
        if (text_.size() > max_record) reject("equilibrium record exceeds 4 MiB");
    }
    // Bound atom allocation while reading untrusted serialized records.
    std::string next() {
        const auto colon = text_.find(':', position_);
        if (colon == std::string_view::npos) reject("missing netstring length");
        std::size_t length = 0;
        const auto size = text_.substr(position_, colon - position_);
        const auto parsed = std::from_chars(size.data(), size.data() + size.size(), length);
        if (size.empty() || (size.size() > 1 && size.front() == '0')
            || parsed.ec != std::errc{} || parsed.ptr != size.data() + size.size()
            || length > 4096 || length >= text_.size() - colon - 1
            || text_[colon + 1 + length] != ',') reject("malformed netstring");
        position_ = colon + length + 2;
        return std::string(text_.substr(colon + 1, length));
    }
    void tag(const char* name) {
        if (next() != name) reject(std::string("unexpected record field: ") + name);
    }
    void text(std::string& value) { value = next(); }
    void number(double& value) {  // a hexadecimal float as the writer prints it
        const auto text = next();
        if (text.empty()) reject("missing record number");
        char* end = nullptr;
        errno = 0;
        value = std::strtod(text.c_str(), &end);
        if (std::any_of(text.begin(), text.end(), [](unsigned char c) { return std::isspace(c) != 0; })
            || end != text.c_str() + text.size() || errno == ERANGE || !std::isfinite(value))
            reject("invalid record number");
    }
    // Bound default component-index and count fields during record admission.
    template<class T> void whole(T& value, std::size_t maximum = 64) {
        const auto text = next();
        std::size_t result = 0;
        const auto parsed = std::from_chars(text.data(), text.data() + text.size(), result);
        if (text.empty() || parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size() || result > maximum)
            reject("invalid record integer");
        value = static_cast<T>(result);
    }
    void count(std::size_t& n, std::size_t maximum) { whole(n, maximum); }
    bool present(const char* name, bool) {  // a named block the writer may omit
        const auto start = position_;
        if (next() == name) return true;
        position_ = start;
        return false;
    }
    void finish() const {
        if (position_ != text_.size()) reject("trailing equilibrium record data");
    }
private:
    // Bound record scanning and memory consumption at the ASL boundary.
    static constexpr std::size_t max_record = 4 * 1024 * 1024;
    std::string_view text_;
    std::size_t position_ = 0;
};

struct EquilibriumContext final {
    std::string record;
    // The record's model at theta0, which owns anchor_result, and the model at another theta.
    std::unique_ptr<HelmholtzModel> model, shifted;
    eq::EquilibriumProblem anchor;
    eq::ContinuationOptions options;
    eq::AslFeedMap map;
    eq::EquilibriumResult anchor_result;
    std::optional<eq::CompiledSystem> compiled;
    eq::EquilibriumResult central;
    std::shared_ptr<eq::SolvedStateActionWorkspace> action_workspace;
    std::vector<double> feed;  // the q of `compiled` and `central`; a null `compiled` is no valid solved state
    struct Rows final {
        std::vector<double> values;
        std::vector<std::vector<double>> gradients, hessians;
        std::vector<char> filled;  // at `feed`: [0] the values, then each gradient column, then each packed Hessian pair
    } rows[2];  // 0: the F C phase-component amounts, then the S solid amounts; 1: the system enthalpy
};

std::string conditioning_warnings(const std::vector<epcsaft::core::NumericalRank>& records) {
    std::ostringstream out; out.precision(17);
    for (const auto& record : records) if (record.first_conditioning_warning) {
        const auto& w = *record.first_conditioning_warning;
        out << "; conditioning warning: " << w.operation << "; basis " << w.basis << "; norm " << w.norm
            << "; method " << w.method << "; cause " << w.cause << "; rank " << w.rank << '/' << w.dimension
            << "; cutoff " << w.threshold << "; condition " << w.condition << "; error " << w.error
            << "; criterion " << w.criterion << ' ' << w.criterion_value << "; bits " << w.precision_bits
            << "; attempt/stage/iteration=" << epcsaft::core::optional_number(w.attempt) << '/' << epcsaft::core::optional_number(w.stage) << '/' << w.iteration << "; point " << w.point_kind;
        for (const auto& [name, values] : {std::pair{"coordinates", &w.coordinates}, {"controls", &w.controls},
            {"row_scale", &w.row_scale}, {"column_scale", &w.column_scale}}) {
            out << "; " << name << '='; for (double value : *values) out << value << ',';
        }
    }
    return out.str();
}

// The refusal of `result` as the context's central state, if any.
std::optional<EquilibriumFailure> central_failure(const EquilibriumContext& context, const eq::EquilibriumResult& result,
                                                 const std::string& observable) {
    const auto failure = [&](ActionStatus status) { return EquilibriumFailure(status, observable, result.message + conditioning_warnings(result.rank_records)
        + (result.first_numerical_failure ? "; first refusal: " + epcsaft::core::describe(*result.first_numerical_failure) : "")); };
    const auto count = context.map.offset.size(), identity = context.map.identity_component;
    if (result.topology_event == eq::TopologyEvent::PhasesCoalesced)
        return failure(ActionStatus::NonuniqueCoalescedPhases);
    if (result.topology_event == eq::TopologyEvent::PhaseDisappeared)
        return failure(ActionStatus::PhaseBoundaryNotDifferentiable);
    if (!result.has_state()) return failure(ActionStatus::CentralSolveUnavailable);
    if (result.topology_event != eq::TopologyEvent::None || result.phase_count != context.anchor.phases.size()
        || result.mole_fractions.size() != result.phase_count * count || result.amounts.size() != result.phase_count
        || result.solid_amounts.size() != context.anchor.solids.size()
        || (result.phase_count == 2 && !(result.mole_fractions[identity] - result.mole_fractions[count + identity] > context.map.minimum_gap)))
        return failure(ActionStatus::PhaseBoundaryNotDifferentiable);
    for (std::size_t s = 0; s < result.solid_amounts.size(); ++s)  // native classification sets an absent amount to 0
        if (!(result.solid_amounts[s] > 0.0))
            return EquilibriumFailure(ActionStatus::PhaseBoundaryNotDifferentiable, observable,
                                     "solid " + context.anchor.solids[s].name + " is absent at the solved state");
    return std::nullopt;
}

std::unique_ptr<EquilibriumContext> prepare_equilibrium(const char* raw, const std::string& observable) {
    RecordReader reader(raw);
    epcsaft::ParameterDocument document;
    std::vector<epcsaft::ActiveParameter> active;
    std::optional<epcsaft::ThermochemistryRecord> thermochemistry;
    auto context = std::make_unique<EquilibriumContext>();
    auto& problem = context->anchor;
    eq::equilibrium_record(reader, document, active, thermochemistry, problem, context->options, context->map);
    reader.finish();
    context->model = std::make_unique<HelmholtzModel>(epcsaft::build_parameters(document), active, thermochemistry, problem.solver_options);
    context->map.validate(problem, context->model->parameters());
    const std::size_t counts[2] = {problem.phases.size() * document.component_ids.size() + problem.solids.size(), 1};
    const std::size_t n = 2 + context->map.anchor.size() + active.size();
    for (int kind = 0; kind < 2; ++kind) {
        auto& rows = context->rows[kind];
        rows.values.resize(counts[kind]);
        rows.gradients.assign(counts[kind], std::vector<double>(n));
        rows.hessians.assign(counts[kind], std::vector<double>(n * (n + 1) / 2));
        rows.filled.assign(1 + n + n * (n + 1) / 2, 0);
    }
    context->anchor_result = eq::solve_equilibrium(*context->model, problem);
    if (auto failure = central_failure(*context, context->anchor_result, observable)) throw *failure;
    context->record = raw;
    return context;
}

ActionStatus continuation_status(eq::ContinuationStop stop) {
    using Stop = eq::ContinuationStop;
    switch (stop) {
    case Stop::invalid_path: case Stop::infeasible: return ActionStatus::DefiningEquationsNotSatisfied;
    case Stop::singular_jacobian: return ActionStatus::RankDeficient;
    case Stop::branch_lost: case Stop::phase_disappeared: case Stop::turning_point: return ActionStatus::PhaseBoundaryNotDifferentiable;
    case Stop::phases_coalesced: return ActionStatus::NonuniqueCoalescedPhases;
    case Stop::left_domain: return ActionStatus::OutsideInterval;
    case Stop::nonfinite: return ActionStatus::Nonfinite;
    default: return ActionStatus::CentralSolveUnavailable;
    }
}

// One native continuation leg from the current central state to `to`; its stop is logged and named after the direct failure.
eq::EquilibriumResult continuation_leg(const EquilibriumContext& context, const HelmholtzModel& model, const eq::EquilibriumProblem& to,
    const std::vector<double>& targets, const std::string& direct, const char* name, const std::string& observable) {
    const auto cause = direct + "; " + name + ": ";
    eq::ContinuationResult result;
    try { result = eq::continue_equilibrium(model, {context.anchor, to, context.central, context.options, targets}); }
    catch (const std::exception& later) { throw EquilibriumFailure(ActionStatus::CentralSolveUnavailable, observable, cause + later.what()); }
    if (!result.success || result.stop != eq::ContinuationStop::reached_target || result.states.empty()) {
        std::fprintf(stderr, "epcsaft_equilibrium diagnostic: %s stop=%s reached=%.17g message=%s\n", name,
                     eq::continuation_stop_name(result.stop), result.reached_path, result.message.c_str());
        if (result.cause_stop) std::fprintf(stderr, "cause_stop=%s\n", eq::continuation_stop_name(*result.cause_stop));
        if (result.cause_message) std::fprintf(stderr, "cause=%s\n", result.cause_message->c_str());
        if (result.event_phase) std::fprintf(stderr, "event_phase=%zu\n", *result.event_phase);
        if (result.event_pair) std::fprintf(stderr, "event_pair=%zu,%zu\n", result.event_pair->first, result.event_pair->second);
        if (result.last_result) std::fprintf(stderr, "solver=%s raw_residual=%.17g classified_residual=%.17g\n",
            result.last_result->solver_status.c_str(), result.last_result->residual_norm(), result.last_result->convergence.residual);
        throw EquilibriumFailure(continuation_status(result.stop), observable, cause + result.message);  // the message leads with its stop
    }
    return result.states.back();
}

// q = (T, P, feed controls, theta). A new q first solves directly in the model at theta; on failure a parameter leg at the
// anchor (T, P, feed) precedes one leg moving T, P and the feed together, since a parameter path moves nothing else. With
// solids there is no continuation: the direct failure keeps its status, and a T outside a solid's law is OutsideInterval.
void select_state(EquilibriumContext& context, const std::vector<double>& q, const std::string& observable) {
    const auto controls = context.map.anchor.size();
    const std::vector<double> feed(q.begin() + 2, q.begin() + 2 + controls), theta(q.begin() + 2 + controls, q.end());
    context.action_workspace.reset();
    context.compiled.reset();  // it points at its model
    if (bitwise_equal(theta, context.model->active_parameter_values())) context.shifted.reset();
    if (!bitwise_equal(theta, (context.shifted ? *context.shifted : *context.model).active_parameter_values())) try {
        context.shifted = model_at(*context.model, theta);  // rebuilt against the current model's own values
    } catch (const std::invalid_argument& error) { throw EquilibriumFailure(ActionStatus::OutsideInterval, observable, error.what()); }
    const auto& model = context.shifted ? *context.shifted : *context.model;
    auto target = context.anchor;
    target.temperature = eq::FreeValue::fixed(q[0]), target.pressure = eq::FreeValue::fixed(q[1]);
    target.feed.values = context.map.amounts(feed);
    for (const double amount : target.feed.values)
        if (!std::isfinite(amount) || amount < 0.0) throw EquilibriumFailure(ActionStatus::OutsideInterval, observable);
    const bool moved = !bitwise_equal(feed, context.map.anchor) || q[0] != context.anchor.temperature.value || q[1] != context.anchor.pressure.value;
    std::optional<eq::CompiledSystem> direct_compiled;
    if (!moved && !context.shifted) context.central = context.anchor_result;
    else {
        std::optional<EquilibriumFailure> refused;
        std::string direct;  // the direct solve's refusal; empty when it succeeds
        try {
            if (target.solids.empty()) context.central = eq::solve_equilibrium_retaining_compiled(model, target, direct_compiled);
            else context.central = eq::solve_equilibrium(model, target);
        }
        catch (const std::exception& error) { direct = error.what(); }
        if (direct.empty() && (refused = central_failure(context, context.central, observable))) direct = refused->what();
        if (!direct.empty()) {
            direct_compiled.reset();
            for (const auto& solid : context.anchor.solids)
                if (q[0] < solid.correlation->temperature_min || q[0] > solid.correlation->temperature_max)
                    throw EquilibriumFailure(ActionStatus::OutsideInterval, observable, "solid " + solid.name + ": " + direct);
            if (!context.anchor.solids.empty()) {
                if (refused) throw *refused;
                throw std::runtime_error(direct);
            }
            const auto leg = [&](const auto& path, const auto& to, const auto& targets, const char* name) {
                context.central = continuation_leg(context, path, to, targets, direct, name, observable);
            };
            context.central = context.anchor_result;
            if (context.shifted) leg(*context.model, context.anchor, theta, "parameter continuation");
            if (moved) leg(model, target, std::vector<double>{}, "continuation");
        }
    }
    if (auto failure = central_failure(context, context.central, observable)) throw *failure;
    context.compiled.emplace(direct_compiled ? std::move(*direct_compiled) : eq::compile_problem(model, target));
    context.feed = q;
    for (auto& rows : context.rows) std::fill(rows.filled.begin(), rows.filled.end(), 0);
}

// One native batch of one kind's rows; an unavailable row fails only its own kind's callback.
// ponytail: a caller reading both kinds at a new q pays both batches; the upgrade is
// one combined batch with per-row failure orders, as PropertyContext::Branch::failures keeps, which could also fill the
// amount rows from an enthalpy call when they are all available.
eq::SolvedStateActionBatch actions(EquilibriumContext& context, const std::vector<int>& coordinates, int kind) {
    eq::SolvedStateActionRequest request;
    if (kind == 1) request.observables.push_back({eq::SolvedStateObservableKind::SystemEnthalpy});
    for (std::size_t phase = 0; kind == 0 && phase < context.anchor.phases.size(); ++phase)
        for (std::size_t component = 0; component < context.map.offset.size(); ++component)
            request.observables.push_back({eq::SolvedStateObservableKind::PhaseComponentAmount, phase, component});
    for (std::size_t s = 0; kind == 0 && s < context.anchor.solids.size(); ++s)
        request.observables.push_back({eq::SolvedStateObservableKind::SolidAmount, s});
    const auto controls = context.map.anchor.size();
    for (const int coordinate : coordinates) {  // a q column (T, P, feed controls, theta); -1 is no direction
        std::vector<double> e(2 + controls + context.model->active_parameters().size(), 0.0);
        if (coordinate >= 0) e.at(coordinate) = 1.0;
        const auto theta = e.begin() + 2 + controls;
        request.directions.push_back({e[0], e[1], context.map.amounts({e.begin() + 2, theta}, true), {theta, e.end()}});
    }
    auto batch = eq::solved_state_actions_prepared(*context.compiled, context.central, request, context.action_workspace);
    for (std::size_t row = 0; row < batch.results.size(); ++row) {
        const auto& item = batch.results[row];
        auto status = item.status;
        if (status == ActionStatus::Available && (!item.value || !item.action
            || !std::isfinite(*item.value) || !std::isfinite(*item.action))) status = ActionStatus::Nonfinite;
        const bool floor_value = coordinates.front() == -1 && status == ActionStatus::DefiningEquationsNotSatisfied
            && !context.central.requested_tolerance_met.value_or(false) && item.value && std::isfinite(*item.value);
        if (status != ActionStatus::Available && !floor_value) {
            const auto& d = batch.diagnostics;
            throw EquilibriumFailure(status, eq::solved_state_observable_id(item.observable), item.diagnostic_message.value_or("")
                + "; power-of-two column-then-row equilibrated defining Jacobian rank "
                + (d.rank ? std::to_string(*d.rank) : "absent") + "/" + std::to_string(d.state_dimension)
                + "; condition " + epcsaft::core::optional_number(d.condition_estimate)
                + "; backward error " + epcsaft::core::optional_number(d.maximum_linear_residual)
                + (item.property_diagnostics && item.property_diagnostics->first_numerical_failure
                    ? "; " + epcsaft::core::describe(*item.property_diagnostics->first_numerical_failure) : ""));
        }
        context.rows[kind].values[row] = *item.value;
    }
    context.rows[kind].filled[0] = 1;
    return batch;
}

// A solver only sees an evaluation error; name the typed cause and the trial real arguments.
double report_failure(arglist* arguments, const char* callback) {
    std::fprintf(stderr, "%s failure: %s controls=", callback, callback_error.c_str());
    for (int i = 0; arguments != nullptr && arguments->ra != nullptr && i < arguments->nr; ++i)
        std::fprintf(stderr, "%s%.17g", i != 0 ? "," : "", arguments->ra[i]);
    std::fputc('\n', stderr);
    invalidate(arguments);
    if (arguments != nullptr) {
        if (arguments->hes != nullptr) callback_error.insert(0, "\"");
        else if (arguments->derivs != nullptr) callback_error.insert(0, "'");
        arguments->Errmsg = callback_error.data();
    }
    return invalid_value;
}

// kind 0: (q, record, phase, component), mol; kind 1: (q, record), the system enthalpy in J; kind 2: (q, record, solid), mol
// of the solid's formula unit. Kinds 0 and 2 read one batch; every kind shares one solved state per record and q. Without
// `second`, no Hessian pair is computed and its read entries keep the NaN written on entry.
double equilibrium_observable(arglist* arguments, int kind, bool second = true) {
    static Prepared<EquilibriumContext> prepared;
    std::lock_guard<std::mutex> lock(prepared.mutex);
    std::string observable = "request";
    if (arguments != nullptr) arguments->Errmsg = nullptr;
    invalidate(arguments);
    try {
        using Kind = eq::SolvedStateObservableKind;
        const Kind kinds[3] = {Kind::PhaseComponentAmount, Kind::SystemEnthalpy, Kind::SolidAmount};
        const int strings[3] = {3, 1, 2};  // the record, then the selectors
        const int n = argument_count(arguments, strings[kind], "equilibrium callback requires real q and the record");
        const auto p = kind != 1 ? selector(arguments, n + 1, kind == 0 ? 2 : 64) : 0, c = kind == 0 ? selector(arguments, n + 2, 64) : 0;
        observable = eq::solved_state_observable_id({kinds[kind], p, c});
        const char* record = string_argument(arguments, n);
        auto& context = prepared.find(record, "equilibrium", [&] { return prepare_equilibrium(record, observable); });
        const auto controls = 2 + context.map.anchor.size();  // T, P and the feed controls
        const auto phases = context.anchor.phases.size(), count = context.map.offset.size();
        if (std::size_t(n) != controls + context.model->active_parameters().size() || c >= count
            || (kind == 0 && p >= phases) || (kind == 2 && p >= context.anchor.solids.size()))
            reject("control/observable dimension mismatch");
        const std::vector<double> q(arguments->ra, arguments->ra + n);  // a parameter may be negative (k_ij)
        for (std::size_t i = 0; i < q.size(); ++i) {
            if (!std::isfinite(q[i])) throw EquilibriumFailure(ActionStatus::Nonfinite, observable);
            if (i < controls && (i < 2 ? q[i] <= 0.0 : q[i] < 0.0)) throw EquilibriumFailure(ActionStatus::OutsideInterval, observable);
        }
        if (!context.compiled || !bitwise_equal(context.feed, q)) select_state(context, q, observable);
        auto& rows = context.rows[kind == 1];  // the amount batch also holds the solid amounts
        const auto row = kind == 2 ? phases * count + p : p * count + c;
        // Each live gradient column and Hessian pair is one batch, computed once per state and kind.
        const auto fill = [&](const std::vector<int>& coordinates, char& filled, std::vector<std::vector<double>>& out, std::size_t entry, real* to) {
            if (std::any_of(coordinates.begin(), coordinates.end(), [&](int i) { return !live(arguments, i); })) return;
            if (!filled) {
                const auto batch = actions(context, coordinates, kind == 1);
                for (std::size_t r = 0; r < rows.values.size(); ++r) out[r][entry] = *batch.results[r].action;
                filled = 1;
            }
            to[entry] = out[row][entry];
        };
        for (int i = 0; arguments->derivs != nullptr && i < n; ++i) fill({i}, rows.filled[1 + i], rows.gradients, i, arguments->derivs);
        for (int j = 0, packed = 0; second && arguments->hes != nullptr && j < n; ++j) for (int i = 0; i <= j; ++i, ++packed)
            fill({i, j}, rows.filled[1 + n + packed], rows.hessians, packed, arguments->hes);
        if (!rows.filled[0]) actions(context, {-1}, kind == 1);
        return rows.values[row];  // On the caller's declared absolute batch basis.
    } catch (const EquilibriumFailure& error) {
        callback_error = error.what();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "epcsaft_equilibrium diagnostic: %s\n", error.what());
        callback_error = "epcsaft_equilibrium:untyped:" + observable + ": " + error.what();
    } catch (...) {
        callback_error = "epcsaft_equilibrium:untyped:" + observable;
    }
    for (const auto& item : prepared.contexts) {  // every record re-solves at its next call
        item->action_workspace.reset();
        item->compiled.reset();
    }
    return report_failure(arguments, "epcsaft_equilibrium");
}

// Pressure-chart property callback (#73): q = (T, P, w_1..w_C, theta_1..theta_m), u = Pw with P = I - zz'/z'z the charge-null projector
// (identity for a neutral record), s = sum(u), x = u/s. Native pressure actions along the O(1), sum- and charge-free
// directions b_j = Pe_j - x sigma_j, sigma_j = 1'Pe_j, supply the value and actions; only dx/dw_j = b_j/s and the exact
// chain d2x/dw_j dw_k = -(sigma_k dx/dw_j + sigma_j dx/dw_k)/s are applied here. Rows: 0 density, 1 total enthalpy,
// 2 + i ln(phi_i), 2 + C + s the solid row l_s. A theta column is one native parameter direction at fixed (T, P, x); x does
// not depend on theta.
struct PropertyContext final {
    std::string record;
    bool require_phase_identity = false;
    std::unique_ptr<HelmholtzModel> model;
    struct Branch final {
        std::vector<double> q, values;
        std::vector<std::vector<double>> gradients, hessians;
        std::vector<std::pair<int, std::string>> failures;  // per row: lowest failing order (3: none) and cause
        int order = -1;
        std::vector<bool> live;  // per q column: whether its derivatives are read
    } solid_rows;  // solid rows need no density root
    struct State final {
        std::size_t phase;
        Branch rows;
        std::optional<epcsaft::DensityRoot> cold_root;
    };
    std::list<State> states;
    std::vector<eq::SolidSpec> solids;
    std::optional<epcsaft::NeutralReferenceDeclaration> reference;
};

std::unique_ptr<PropertyContext> prepare_property(const char* raw) {
    RecordReader reader(raw);
    epcsaft::ParameterDocument document;
    std::vector<epcsaft::ActiveParameter> active;
    std::optional<epcsaft::ThermochemistryRecord> thermochemistry;
    epcsaft::policy::SolverOptions solver_options;
    auto context = std::make_unique<PropertyContext>();
    eq::property_record(reader, context->require_phase_identity, document, active, thermochemistry, solver_options,
                        context->solids, context->reference);
    reader.finish();
    context->model = std::make_unique<HelmholtzModel>(
        epcsaft::build_parameters(document), active, thermochemistry, solver_options);
    eq::admit_solid_declarations(*context->model, context->reference, context->solids);
    context->record = raw;
    return context;
}

const PropertyContext::Branch& evaluate_property(PropertyContext& context, std::size_t phase, const std::vector<double>& q, const std::vector<bool>& live, int order) {
    auto& states = context.states;
    const auto found = std::find_if(states.begin(), states.end(), [&](const auto& state) {
        return state.phase == phase && state.rows.live == live && bitwise_equal(state.rows.q, q);
    });
    if (found != states.end() && found->rows.order >= order) {
        states.splice(states.begin(), states, found);
        return states.front().rows;
    }
    auto state = found == states.end() ? PropertyContext::State{phase, {}, {}} : *found;
    auto& branch = state.rows;
    using Status = epcsaft::PropertyActionStatus;
    const std::size_t c = context.model->parameters().component_count(), rows = 2 + c, n = q.size();
    const auto active = context.model->active_parameters();
    const std::vector<double> theta(q.begin() + rows, q.end());
    if (!bitwise_equal(theta, context.model->active_parameter_values())) context.model = model_at(*context.model, theta);
    const auto& z = context.model->parameters().charge_number;
    const double zz = std::inner_product(z.begin(), z.end(), z.begin(), 0.0);
    std::vector<double> u(q.begin() + 2, q.begin() + 2 + c), x(c), sigma(c, 1.0);
    if (zz > 0.0) {  // admit within the native tolerance, then remove only the rounding below it
        const double total = std::accumulate(u.begin(), u.end(), 0.0);
        for (std::size_t i = 0; i < c; ++i) x[i] = u[i] / total;
        context.model->validate_electroneutrality(x);
        const double zw = std::inner_product(z.begin(), z.end(), u.begin(), 0.0), z1 = std::accumulate(z.begin(), z.end(), 0.0);
        for (std::size_t i = 0; i < c; ++i) u[i] -= z[i] * zw / zz, sigma[i] = 1.0 - z[i] * z1 / zz;
    }
    const double sum = std::accumulate(u.begin(), u.end(), 0.0);
    for (std::size_t i = 0; i < c; ++i) x[i] = u[i] / sum;
    std::vector<epcsaft::PressureActionDirection> basis(rows, {0.0, 0.0, std::vector<double>(c, 0.0), {}});
    basis[0].temperature = 1.0;
    basis[1].pressure = 1.0;
    for (std::size_t j = 0; j < c; ++j)
        for (std::size_t i = 0; i < c; ++i)
            basis[2 + j].composition[i] = (i == j ? 1.0 : 0.0) - (zz > 0.0 ? z[i] * z[j] / zz : 0.0) - x[i] * sigma[j];
    for (const auto& a : active) basis.push_back({0.0, 0.0, std::vector<double>(c, 0.0), {{a.family, a.indices, 1.0}}});
    const std::vector<epcsaft::PropertyObservable> observables = {
        epcsaft::PropertyObservable::TotalEnthalpy, epcsaft::PropertyObservable::LogFugacityCoefficient};
    const auto cause = [](Status status, std::string fallback) {
        if (status == Status::PressureSingular) return std::string("nonmechanical_root");
        return status == Status::SupportChange || status == Status::IonicInsertion ? "support_insertion" : fallback;
    };
    // The first, exhaustive root search at this q owns branch identity and every value; later directions start from its root.
    std::optional<double> anchor;
    if (branch.order < 0) {
        branch = {q, std::vector<double>(rows), std::vector<std::vector<double>>(rows, std::vector<double>(n)),
            std::vector<std::vector<double>>(rows, std::vector<double>(n * (n + 1) / 2)),
            std::vector<std::pair<int, std::string>>(rows, {3, ""}), -1, live};
    } else if (branch.order != 0 || !live[0]) anchor = branch.values[0];
    const auto requested = phase == 0 ? epcsaft::DensityBranch::Liquid : epcsaft::DensityBranch::Vapor;
    // With live T, its first action uses the cold root, even after a value-only request.
    context.model->pressure_actions(q[0], q[1], x, requested, observables, anchor,
        !anchor && live[0] && (order == 0 || state.cold_root) ? &state.cold_root : nullptr,
        [&](const HelmholtzModel::PressureAction& evaluate) {
            const auto run = [&](const std::vector<epcsaft::PressureActionDirection>& directions) {
                const auto batch = evaluate(directions);
                if (batch.density_status != Status::Available || !batch.molar_density) {
                    const auto& m = batch.density_message;
                    std::string root = batch.density_convergence ? "inner_iteration_failure" : "unavailable_derivative";
                    if (m.find("ambiguous") != std::string::npos) root = "ambiguous_branch";
                    if (m.find("no stable density root") != std::string::npos || m.find("density domain is empty") != std::string::npos)
                        root = "no_density_root";
                    throw std::runtime_error("epcsaft_property:" + cause(batch.density_status, root) + ":"
                        + std::to_string(static_cast<int>(batch.density_status)) + ":density:"
                        + (batch.diagnostics.first_numerical_failure ? epcsaft::core::describe(*batch.diagnostics.first_numerical_failure)
                            : batch.density_convergence ? epcsaft::core::describe(*batch.density_convergence) : "") + "; " + m);
                }
                const auto returned = batch.diagnostics.returned_branch;
                if ((!anchor && context.require_phase_identity && returned == epcsaft::DensityBranch::Unspecified)
                    || (returned != epcsaft::DensityBranch::Unspecified && returned != requested))
                    throw std::runtime_error("epcsaft_property:branch_mismatch:requested=" + std::to_string(static_cast<int>(requested))
                        + ":returned=" + std::to_string(static_cast<int>(returned)) + ":T=" + epcsaft::core::optional_number(q[0])
                        + ":P=" + epcsaft::core::optional_number(q[1]) + ":stable_root_count="
                        + (batch.diagnostics.stable_root_count ? std::to_string(*batch.diagnostics.stable_root_count) : "absent"));
                const auto fail = [&](std::size_t row, int k, Status status) {
                    if (k < branch.failures[row].first)
                        branch.failures[row] = {k, "epcsaft_property:" + cause(status, status == Status::MissingThermochemistry ? "MissingThermochemistry" : status == Status::OutsideInterval ? "OutsideInterval" : "unavailable_derivative") + ":"
                            + std::to_string(static_cast<int>(status)) + ":" + std::to_string(row)
                            + ((status == Status::ThermalBoundary || status == Status::OutsideInterval) ? ": thermochemistry interval at T=" + epcsaft::core::optional_number(q[0]) + " K" : "")
                            + ((status == Status::OutsideInterval || status == Status::MissingThermochemistry) ? "" : "; site rank " + std::to_string(batch.diagnostics.site_rank)
                            + "; pivot ratio " + epcsaft::core::optional_number(batch.diagnostics.site_pivot_ratio)
                            + "; association residual " + epcsaft::core::optional_number(batch.diagnostics.association_residual))
                            + (batch.diagnostics.first_numerical_failure ? "; " + epcsaft::core::describe(*batch.diagnostics.first_numerical_failure) : "")};
                };
                std::vector<double> actions(rows);
                const bool first = !anchor && branch.order < 0;
                if (!anchor) anchor = *batch.molar_density;
                if (first) branch.values[0] = *anchor;
                if (!directions.empty()) actions[0] = directions.size() == 1 ? batch.density_action.at(0) : batch.density_mixed_action.value();
                for (std::size_t row = 1; row < rows; ++row) {
                    const auto& item = batch.results[row - 1];
                    if (first && item.value) branch.values[row] = *item.value;  // an owned thermochemistry endpoint keeps its value
                    if (first && !item.value) fail(row, 0, item.status);
                    actions[row] = item.action.value_or(0.0);
                    if (!directions.empty() && (item.status != Status::Available || !item.action)) fail(row, static_cast<int>(directions.size()), item.status);
                }
                return actions;
            };
            for (std::size_t k = 0; branch.order < 1 && k < (order >= 1 ? n : 1); ++k) {
                if (k == 0 && (order == 0 || !live[0])) {
                    if (!anchor) run({});
                    continue;
                }
                if (!live[k]) continue;
                const auto actions = run({basis[k]});
                for (std::size_t row = 0; row < rows; ++row) branch.gradients[row][k] = actions[row] / std::pow(sum, k >= 2 && k < rows);
            }
            std::size_t packed = 0;
            for (std::size_t j = 0; order >= 2 && j < n; ++j)
                for (std::size_t i = 0; i <= j; ++i, ++packed) {
                    if (!live[i] || !live[j]) continue;
                    const auto actions = run({basis[i], basis[j]});
                    for (std::size_t row = 0; row < rows; ++row)
                        branch.hessians[row][packed] = actions[row] / std::pow(sum, (i >= 2 && i < rows) + (j >= 2 && j < rows))
                            - (i >= 2 && j < rows ? (sigma[j - 2] * branch.gradients[row][i] + sigma[i - 2] * branch.gradients[row][j]) / sum : 0.0);
                }
        });
    branch.order = order;
    if (found != states.end()) states.erase(found);
    states.push_front(std::move(state));
    if (states.size() > 128) states.pop_back();
    return states.front().rows;
}

// EqID solid_standard_gibbs on the record basis: l_s = ln K^rho_s(T, P0) + V_s (P - P0)/RT - dnu ln(P/(RT rho0)), rho0 =
// 1 mol/m^3, so sum_i nu_i (ln x_i + ln phi_i) = l_s at saturation. It depends on (T, P, theta) only.
void evaluate_solids(PropertyContext& context, const std::vector<double>& q, const std::vector<bool>& live, int order) {
    auto& rows = context.solid_rows;
    if (bitwise_equal(rows.q, q) && rows.live == live && rows.order >= order) return;
    const std::size_t S = context.solids.size(), n = q.size(), first = n - context.model->active_parameters().size();
    const std::vector<double> theta(q.begin() + first, q.end());
    if (!bitwise_equal(theta, context.model->active_parameter_values())) context.model = model_at(*context.model, theta);
    std::vector<std::size_t> parameters;  // the live parameters: their first-order jets, and pairs at order 2
    for (std::size_t k = first; order >= 1 && k < n; ++k) if (live[k]) parameters.push_back(k - first);
    rows = {q, std::vector<double>(S), std::vector<std::vector<double>>(S, std::vector<double>(n)),
            std::vector<std::vector<double>>(S, std::vector<double>(n * (n + 1) / 2)), std::vector<std::pair<int, std::string>>(S, {3, ""}), order};
    rows.live = live;
    const double T = q[0], P = q[1], RT = epcsaft::GAS_CONSTANT_J_PER_MOL_K * T;
    eq::ReferenceJets jets;
    jets.temperature_order = live[0] ? std::max(order, 0) : 0;
    jets.parameters = parameters;
    if (order >= 2) jets.pairs = parameters;
    for (std::size_t s = 0; s < S; ++s) {
        const auto& solid = context.solids[s];
        const auto row = std::to_string(first + s) + ": solid " + solid.name + ": ";
        const auto fail = [&](int k, const char* cause, const std::string& detail) {
            if (k < rows.failures[s].first) rows.failures[s] = {k, "epcsaft_property:" + (cause + (":" + row)) + detail};
        };
        const double P0 = *solid.reference->reference_pressure_pa;
        eq::SolidLogK law;
        try {  // the reference is evaluated at its fixed source pressure, so it is independent of P
            law = eq::solid_log_k(*context.model, T, P0, context.reference, {solid}, jets)[0];
        } catch (const std::exception& error) {
            const bool outside = T < solid.correlation->temperature_min || T > solid.correlation->temperature_max;
            fail(0, outside ? "OutsideInterval" : "unavailable_derivative", error.what());
            continue;
        }
        const auto& C = law.reference;
        const double V = solid.molar_volume_m3_per_mol, dnu = std::accumulate(solid.stoichiometry.begin(), solid.stoichiometry.end(), 0.0);
        rows.values[s] = law.ln_k + law.scale + *C.terminal_contraction + V * (P - P0) / RT - dnu * std::log(P / RT);
        // A requested column the reference owner did not admit fails the row at its order.
        const auto column = [&](int k, std::size_t a, std::optional<std::size_t> i = {}, std::optional<std::size_t> j = {}) {
            const auto v = C.column(a, 0, i, j);
            if (v.state == epcsaft::ReferenceColumnState::Available) return v.value;
            fail(k, "unavailable_derivative", "reference column d^" + std::to_string(a) + "/dT; " + C.message);
            return 0.0;
        };
        auto& g = rows.gradients[s];
        const auto h = [&](std::size_t i, std::size_t j) -> double& { return rows.hessians[s][j * (j + 1) / 2 + i]; };  // i <= j
        if (order >= 1) {
            if (live[0]) g[0] = *law.slope + column(1, 1) - V * (P - P0) / (RT * T) + dnu / T;
            if (live[1]) g[1] = V / RT - dnu / P;
            for (std::size_t k = first; k < n; ++k) if (live[k]) g[k] = column(1, 0, k - first);
        }
        if (order < 2) continue;
        if (live[0]) h(0, 0) = *law.curvature + column(2, 2) + 2.0 * V * (P - P0) / (RT * T * T) - dnu / (T * T);
        if (live[0] && live[1]) h(0, 1) = -V / (RT * T);
        if (live[1]) h(1, 1) = dnu / (P * P);
        for (std::size_t k = first; k < n; ++k) {
            if (live[0] && live[k]) h(0, k) = column(2, 1, k - first);
            for (std::size_t l = first; l <= k; ++l) if (live[l] && live[k]) h(l, k) = column(2, 0, l - first, k - first);
        }
    }
}

double pressure_property(arglist* arguments, bool second) {
    static Prepared<PropertyContext> prepared;
    std::lock_guard<std::mutex> lock(prepared.mutex);
    if (arguments != nullptr) arguments->Errmsg = nullptr;
    invalidate(arguments);
    try {
        const int n = argument_count(arguments, 3, "property callback requires real arguments and three strings");
        const char* record = string_argument(arguments, n);
        auto& context = prepared.find(record, "property", [&] { return prepare_property(record); });
        const std::size_t c = context.model->parameters().component_count();
        if (std::size_t(n) != 2 + c + context.model->active_parameters().size()) reject("argument count differs");
        const auto row = selector(arguments, n + 1, 2 + c + context.solids.size()), phase = selector(arguments, n + 2, 2);
        const std::vector<double> q(arguments->ra, arguments->ra + n);
        // Zero weights are admitted: the component enters at its exact infinite-dilution limit.
        for (std::size_t i = 0; i < q.size(); ++i)  // a parameter may be negative (k_ij)
            if (!std::isfinite(q[i]) || (i < 2 + c && q[i] < 0.0) || (i < 2 && q[i] == 0.0))
                reject("T and P must be finite and positive, weights finite and nonnegative, parameters finite");
        if (!(std::accumulate(q.begin() + 2, q.begin() + 2 + c, 0.0) > 0.0)) reject("weights must have a positive sum");
        const int order = second && arguments->hes != nullptr ? 2 : (arguments->derivs != nullptr ? 1 : 0);
        const bool solid = row >= 2 + c;
        std::vector<bool> used(n);
        for (int i = 0; i < n; ++i) used[i] = live(arguments, i);
        if (solid) evaluate_solids(context, q, used, order);
        const auto& branch = solid ? context.solid_rows : evaluate_property(context, phase, q, used, order);
        const auto r = solid ? row - 2 - c : row;
        if (branch.failures[r].first <= order) throw std::runtime_error(branch.failures[r].second);
        for (int j = 0, packed = 0; j < n; ++j) {  // only the entries ASL reads
            if (arguments->derivs != nullptr && used[j]) arguments->derivs[j] = branch.gradients[r][j];
            for (int i = 0; i <= j; ++i, ++packed)
                if (second && arguments->hes != nullptr && used[i] && used[j]) arguments->hes[packed] = branch.hessians[r][packed];
        }
        return branch.values[r];
    } catch (const std::invalid_argument& error) {
        callback_error = std::string("epcsaft_property:invalid_input:") + error.what();
    } catch (const std::exception& error) {
        const std::string what = error.what();
        callback_error = what.rfind("epcsaft_property:", 0) == 0 ? what : "epcsaft_property:native_failure:" + what;
    }
    return report_failure(arguments, "epcsaft_property");
}

}  // namespace

// A `_first_order` name computes no Hessian pair, for a solver that reads no Hessian; it shares the exact name's solved states.
extern "C" void funcadd_ASL(AmplExports* exports) {
    if (exports == nullptr || exports->Addfunc == nullptr) return;
    const std::tuple<const char*, rfunc, int> functions[] = {
        {"epcsaft_equilibrium_phase_component_amount", [](arglist* a) { return equilibrium_observable(a, 0); }, -5},
        {"epcsaft_equilibrium_system_enthalpy", [](arglist* a) { return equilibrium_observable(a, 1); }, -3},
        {"epcsaft_equilibrium_solid_amount", [](arglist* a) { return equilibrium_observable(a, 2); }, -4},
        {"epcsaft_equilibrium_phase_component_amount_first_order", [](arglist* a) { return equilibrium_observable(a, 0, false); }, -5},
        {"epcsaft_equilibrium_system_enthalpy_first_order", [](arglist* a) { return equilibrium_observable(a, 1, false); }, -3},
        {"epcsaft_equilibrium_solid_amount_first_order", [](arglist* a) { return equilibrium_observable(a, 2, false); }, -4},
        {"epcsaft_pressure_property", [](arglist* a) { return pressure_property(a, true); }, -7},
        {"epcsaft_pressure_property_first_order", [](arglist* a) { return pressure_property(a, false); }, -7}};
    for (const auto& [name, function, count] : functions)
        exports->Addfunc(name, function, FUNCADD_REAL_VALUED | FUNCADD_STRING_ARGS, count, nullptr, exports);
}
