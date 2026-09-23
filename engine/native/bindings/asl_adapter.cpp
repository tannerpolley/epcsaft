#include "properties/helmholtz_model.hpp"
#include "surface/parameter_builder.hpp"
#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium_record.hpp"

#include <charconv>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>

// Use the process stderr for retained diagnostics, including direct Pyomo calls without AE stdio hooks.
#define NO_STDIO1
#define No_AE_redefs
#include "../../build/subprojects/packagefiles/asl/funcadd.h"

namespace {

using epcsaft::DerivativeOrder;
using epcsaft::NeutralProperty;
using epcsaft::NeutralPropertyDerivatives;

struct PreparedModel final {
    std::mutex mutex;
    std::tuple<std::string, std::string, std::string, std::string> inputs;
    std::unique_ptr<epcsaft::HelmholtzModel> model;
};

PreparedModel& prepared_model() {
    static PreparedModel result;
    return result;
}

thread_local std::string callback_error;

double invalid_value() {
    return std::numeric_limits<double>::quiet_NaN();
}

void invalidate(arglist* arguments) {
    if (arguments == nullptr) return;
    if (arguments->derivs != nullptr) {
        for (int index = 0; index < arguments->nr; ++index) {
            arguments->derivs[index] = invalid_value();
        }
    }
    if (arguments->hes != nullptr) {
        const int count = arguments->nr * (arguments->nr + 1) / 2;
        for (int index = 0; index < count; ++index) {
            arguments->hes[index] = invalid_value();
        }
    }
}

[[noreturn]] void reject(const std::string& message) {
    throw std::invalid_argument(message);
}

double parse_decimal(const char* value, const char* name) {
    if (value == nullptr || *value == '\0') reject(std::string("missing ") + name);
    const std::string text(value);
    for (const char character : text) {
        if (std::isspace(static_cast<unsigned char>(character)) != 0) {
            reject(std::string("invalid ") + name);
        }
    }
    std::size_t end = 0;
    double result = 0.0;
    try {
        result = std::stod(text, &end);
    } catch (const std::exception&) {
        reject(std::string("invalid ") + name);
    }
    if (end != text.size() || !std::isfinite(result)) {
        reject(std::string("invalid ") + name);
    }
    return result;
}

void validate_arguments(arglist* arguments) {
    if (arguments == nullptr || arguments->n != 6 || arguments->nr != 2
        || arguments->at == nullptr || arguments->ra == nullptr
        || arguments->sa == nullptr) {
        reject("epcsaft neutral function requires six arguments (two real, four string)");
    }
    for (int position = 0; position < arguments->n; ++position) {
        const int expected = position < 2 ? position : -(position - 1);
        if (arguments->at[position] != expected) reject("invalid neutral argument map");
    }
}

double real_argument(arglist* arguments, int position) {
    const int index = arguments->at[position];
    if (index < 0 || index >= arguments->nr) reject("invalid real argument map");
    return arguments->ra[index];
}

const char* string_argument(arglist* arguments, int position) {
    const int index = arguments->at[position];
    if (index >= 0) reject("invalid string argument map");
    const char* value = arguments->sa[-index - 1];
    if (value == nullptr) reject("missing string argument");
    return value;
}

epcsaft::HelmholtzModel& model_for(
    const std::tuple<std::string, std::string, std::string, std::string>& inputs
) {
    auto& prepared = prepared_model();
    if (prepared.model != nullptr) {
        if (prepared.inputs != inputs) reject("immutable neutral parameter tuple changed");
        return *prepared.model;
    }

    const auto [m, sigma, epsilon, molar_mass] = inputs;
    epcsaft::ParameterDocument document;
    document.component_ids = {"component-0"};
    document.base_choice = "pc-saft";
    document.association_choice = "none";
    document.electrolyte_choice = "none";
    document.ionic_dispersion_choice = "";
    document.permittivity_choice = "none";
    auto add = [&](const char* field, double magnitude, const char* unit) {
        epcsaft::ParameterEntry entry;
        entry.identity = std::string("component/0/") + field;
        entry.field = field;
        entry.index = {0};
        entry.magnitude = magnitude;
        entry.unit = unit;
        document.entries.push_back(std::move(entry));
    };
    add("segment_count", parse_decimal(m.c_str(), "segment_count"), "dimensionless");
    add("sigma_constant", parse_decimal(sigma.c_str(), "sigma_constant"), "angstrom");
    add(
        "dispersion_energy_over_k",
        parse_decimal(epsilon.c_str(), "dispersion_energy_over_k"),
        "kelvin"
    );
    add("molar_mass", parse_decimal(molar_mass.c_str(), "molar_mass"), "kilogram / mole");
    add("charge_number", 0.0, "elementary-charge");
    prepared.model = std::make_unique<epcsaft::HelmholtzModel>(
        epcsaft::build_parameters(document)
    );
    prepared.inputs = inputs;
    return *prepared.model;
}

double evaluate(
    arglist* arguments,
    NeutralProperty property
) {
    validate_arguments(arguments);
    const auto inputs = std::make_tuple(
        std::string(string_argument(arguments, 2)),
        std::string(string_argument(arguments, 3)),
        std::string(string_argument(arguments, 4)),
        std::string(string_argument(arguments, 5))
    );
    const double temperature = real_argument(arguments, 0);
    const double density = real_argument(arguments, 1);
    auto& prepared = prepared_model();
    std::lock_guard<std::mutex> lock(prepared.mutex);
    const auto order = arguments->hes != nullptr
        ? DerivativeOrder::Hessian
        : (arguments->derivs != nullptr ? DerivativeOrder::Gradient : DerivativeOrder::Value);
    const NeutralPropertyDerivatives result = model_for(inputs).neutral_property(
        temperature, density, property, order
    );
    if (arguments->derivs != nullptr) {
        for (int position = 0; position < arguments->n; ++position) {
            const int index = arguments->at[position];
            if (index >= 0) {
                if (index >= arguments->nr) reject("invalid derivative map");
                arguments->derivs[index] = result.gradient.at(position);
            }
        }
    }
    if (arguments->hes != nullptr) {
        arguments->hes[0] = result.hessian.at(0);
        arguments->hes[1] = result.hessian.at(1);
        arguments->hes[2] = result.hessian.at(3);
    }
    return result.value;
}

double callback(arglist* arguments, NeutralProperty property) {
    if (arguments != nullptr) arguments->Errmsg = nullptr;
    invalidate(arguments);
    try {
        return evaluate(arguments, property);
    } catch (const std::exception& error) {
        callback_error = arguments != nullptr && arguments->hes != nullptr
            ? std::string("\"") + error.what()
            : (arguments != nullptr && arguments->derivs != nullptr
                ? std::string("'") + error.what()
                : error.what());
        if (arguments != nullptr) arguments->Errmsg = callback_error.data();
        return invalid_value();
    } catch (...) {
        callback_error = "unknown epcsaft neutral callback error";
        if (arguments != nullptr) arguments->Errmsg = callback_error.data();
        return invalid_value();
    }
}

double residual_helmholtz(arglist* arguments) {
    return callback(arguments, NeutralProperty::ResidualHelmholtz);
}

double pressure(arglist* arguments) {
    return callback(arguments, NeutralProperty::Pressure);
}

double recordings(arglist* arguments) {
    if (arguments == nullptr || arguments->n != 0 || arguments->nr != 0) {
        callback_error = "epcsaft neutral recordings requires zero arguments";
        if (arguments != nullptr) arguments->Errmsg = callback_error.data();
        return invalid_value();
    }
    auto& prepared = prepared_model();
    std::lock_guard<std::mutex> lock(prepared.mutex);
    return prepared.model == nullptr
        ? 0.0
        : static_cast<double>(prepared.model->counters().recordings);
}

namespace eq = epcsaft::equilibrium;
using ActionStatus = eq::SolvedStateActionStatus;

class EquilibriumFailure final : public std::runtime_error {
public:
    EquilibriumFailure(ActionStatus status, const std::string& observable)
        : std::runtime_error(std::string("epcsaft_equilibrium:")
              + eq::solved_state_status_name(status) + ":" + observable) {}
};

// Record admission reconstructs native inputs; the Surface builder owns parameter rules.
class RecordReader final {
public:
    explicit RecordReader(const char* raw) : text_(raw, strnlen(raw, max_record + 1)) {
        if (text_.size() > max_record) reject("equilibrium record exceeds 4 MiB");
    }
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
    void field(const char* name) {
        if (next() != name) reject(std::string("unexpected record field: ") + name);
    }
    std::size_t integer(std::size_t maximum = 64) {
        const auto value = next();
        std::size_t result = 0;
        const auto parsed = std::from_chars(value.data(), value.data() + value.size(), result);
        if (value.empty() || parsed.ec != std::errc{}
            || parsed.ptr != value.data() + value.size() || result > maximum)
            reject("invalid record integer");
        return result;
    }
    double number() { return parse_decimal(next().c_str(), "record number"); }
    template<class Read> auto vector(Read read, std::size_t maximum = 64) {
        std::vector<decltype(read())> result;
        const auto count = integer(maximum);
        result.reserve(count);
        for (std::size_t i = 0; i < count; ++i) result.push_back(read());
        return result;
    }
    void finish() const {
        if (position_ != text_.size()) reject("trailing equilibrium record data");
    }
private:
    static constexpr std::size_t max_record = 4 * 1024 * 1024;
    std::string_view text_;
    std::size_t position_ = 0;
};

struct EquilibriumContext final {
    std::string record;
    std::unique_ptr<epcsaft::HelmholtzModel> model;
    eq::EquilibriumProblem anchor;
    eq::ContinuationOptions options;
    eq::AslFeedMap map;
    eq::EquilibriumResult anchor_result;
    std::unique_ptr<eq::CompiledSystem> compiled;
    eq::EquilibriumResult central;
    std::vector<double> feed, values;
    std::vector<std::vector<double>> gradients, hessians;
    int order = -1;
};

struct PreparedEquilibrium final {
    // ponytail: one immutable model per process; qualify multiple contexts before extending.
    std::mutex mutex;
    std::unique_ptr<EquilibriumContext> context;
};

PreparedEquilibrium& prepared_equilibrium() {
    static PreparedEquilibrium result;
    return result;
}

void require_central(const EquilibriumContext& context, const eq::EquilibriumResult& result, const std::string& observable) {
    const auto count = context.map.offset.size(), identity = context.map.identity_component;
    if (result.topology_event == eq::TopologyEvent::PhasesCoalesced)
        throw EquilibriumFailure(ActionStatus::NonuniqueCoalescedPhases, observable);
    if (result.topology_event == eq::TopologyEvent::PhaseDisappeared)
        throw EquilibriumFailure(ActionStatus::PhaseBoundaryNotDifferentiable, observable);
    if (!result.has_state()) throw EquilibriumFailure(ActionStatus::CentralSolveUnavailable, observable);
    if (!result.requested_tolerance_met.value_or(false))
        throw EquilibriumFailure(ActionStatus::DefiningEquationsNotSatisfied, observable);
    if (result.topology_event != eq::TopologyEvent::None || result.phase_count != 2
        || result.mole_fractions.size() != 2 * count || result.amounts.size() != 2
        || !(result.mole_fractions[identity] - result.mole_fractions[count + identity] > context.map.minimum_gap))
        throw EquilibriumFailure(ActionStatus::PhaseBoundaryNotDifferentiable, observable);
    if (!std::isfinite(result.reduced_hessian_minimum) || result.reduced_hessian_minimum <= 0.0)
        throw EquilibriumFailure(ActionStatus::CentralPropertyUnavailable, observable);
}

std::unique_ptr<EquilibriumContext> prepare_equilibrium(const char* raw, const std::string& observable) {
    RecordReader reader(raw);
    reader.field("epcsaft.asl.equilibrium.v3");
    epcsaft::ParameterDocument document;
    reader.field("component_ids");
    document.component_ids = reader.vector([&] { return reader.next(); });
    reader.field("site_component");
    document.site_component = reader.vector([&] { return static_cast<int>(reader.integer()); });
    reader.field("site_multiplicity");
    document.site_multiplicity = reader.vector([&] { return static_cast<int>(reader.integer()); });
    reader.field("base_choice");
    document.base_choice = reader.next();
    reader.field("association_choice");
    document.association_choice = reader.next();
    reader.field("electrolyte_choice");
    document.electrolyte_choice = reader.next();
    reader.field("ionic_dispersion_choice");
    document.ionic_dispersion_choice = reader.next();
    reader.field("permittivity_choice");
    document.permittivity_choice = reader.next();
    reader.field("parameter_entries");
    document.entries = reader.vector([&] {
        epcsaft::ParameterEntry entry;
        entry.identity = reader.next();
        entry.field = reader.next();
        entry.index = reader.vector([&] { return reader.integer(); });
        entry.magnitude = reader.number();
        entry.unit = reader.next();
        return entry;
    }, 4096);
    auto context = std::make_unique<EquilibriumContext>();
    auto& problem = context->anchor;
    reader.field("temperature");
    problem.temperature = eq::FreeValue::fixed(reader.number());
    reader.field("pressure");
    problem.pressure = eq::FreeValue::fixed(reader.number());
    reader.field("feed_kind");
    if (reader.next() != "amounts") reject("equilibrium callback requires absolute amounts");
    reader.field("feed_amounts");
    problem.feed.values = reader.vector([&] { return reader.number(); });
    reader.field("phases");
    problem.phases = reader.vector([&] {
        eq::PhaseSpec phase;
        phase.name = reader.next();
        if (reader.next() != "liquid") reject("equilibrium callback requires liquid phases");
        phase.support = reader.vector([&] { return reader.integer(); });
        phase.support_declared = reader.integer(1) != 0;
        if (reader.next() != "free") reject("equilibrium callback requires free phase amounts");
        phase.amount.fraction = reader.number();
        phase.composition_guess = reader.vector([&] { return reader.number(); });
        const auto packing = reader.vector([&] { return reader.number(); }, 1);
        if (!packing.empty()) phase.packing_guess = packing.front();
        return phase;
    });
    reader.field("reactions");
    problem.reactions = reader.vector([&] {
        eq::ReactionSpec r;
        r.name = reader.next();
        r.stoichiometry = reader.vector([&] { return reader.number(); });
        const auto target = reader.next();
        if (target == "K") r.equilibrium_constant = reader.number();
        else if (target == "affinity") r.affinity_j_per_mol = reader.number();
        else reject("unsupported reaction target");
        return r;
    });
    reader.field("standard_gibbs_j_per_mol");
    problem.standard_gibbs_j_per_mol = reader.vector([&] { return reader.number(); });
    reader.field("topology_declared");
    problem.topology_declared = reader.integer(1) != 0;
    reader.field("row_tolerance");
    problem.row_tolerance = reader.number();
    reader.field("continuation_options");
    context->options.initial_step = reader.number();
    context->options.minimum_step = reader.number();
    context->options.max_steps = reader.integer();
    context->options.max_retries_per_step = reader.integer();
    reader.field("feed_offset");
    context->map.offset = reader.vector([&] { return reader.number(); });
    reader.field("feed_matrix");
    context->map.matrix = reader.vector([&] { return reader.vector([&] { return reader.number(); }); });
    reader.field("anchor_controls");
    context->map.anchor = reader.vector([&] { return reader.number(); });
    reader.field("phase_identity");
    context->map.identity_component = reader.integer();
    context->map.minimum_gap = reader.number();
    reader.finish();
    context->model = std::make_unique<epcsaft::HelmholtzModel>(epcsaft::build_parameters(document));
    context->map.validate(problem, context->model->parameters());
    const auto rows = 2 * document.component_ids.size(), n = context->map.anchor.size();
    context->values.resize(rows);
    context->gradients.assign(rows, std::vector<double>(n));
    context->hessians.assign(rows, std::vector<double>(n * (n + 1) / 2));
    context->anchor_result = eq::solve_equilibrium(*context->model, problem);
    require_central(*context, context->anchor_result, observable);
    context->record = raw;
    return context;
}

ActionStatus continuation_status(eq::ContinuationStop stop) {
    using Stop = eq::ContinuationStop;
    switch (stop) {
    case Stop::invalid_path: case Stop::infeasible: return ActionStatus::DefiningEquationsNotSatisfied;
    case Stop::singular_jacobian: return ActionStatus::RankDeficient;
    case Stop::branch_lost: case Stop::phase_disappeared: return ActionStatus::PhaseBoundaryNotDifferentiable;
    case Stop::phases_coalesced: return ActionStatus::NonuniqueCoalescedPhases;
    case Stop::left_domain: return ActionStatus::OutsideInterval;
    case Stop::nonfinite: return ActionStatus::Nonfinite;
    default: return ActionStatus::CentralSolveUnavailable;
    }
}

void select_feed(EquilibriumContext& context, const std::vector<double>& feed,
                 const std::string& observable) {
    auto target = context.anchor;
    target.feed.values = context.map.amounts(feed);
    for (const double amount : target.feed.values)
        if (!std::isfinite(amount) || amount < 0.0) throw EquilibriumFailure(ActionStatus::OutsideInterval, observable);
    if (feed == context.map.anchor) context.central = context.anchor_result;
    else {
        const eq::ContinuationRequest request{context.anchor, target, context.anchor_result, context.options};
        const auto result = eq::continue_equilibrium(*context.model, request);
        if (!result.success || result.stop != eq::ContinuationStop::reached_target || result.states.empty()) {
            std::fprintf(stderr, "epcsaft_equilibrium diagnostic: stop=%s reached=%.17g message=%s\n",
                         eq::continuation_stop_name(result.stop), result.reached_path, result.message.c_str());
            if (result.cause_stop) std::fprintf(stderr, "cause_stop=%s\n", eq::continuation_stop_name(*result.cause_stop));
            if (result.cause_message) std::fprintf(stderr, "cause=%s\n", result.cause_message->c_str());
            if (result.event_phase) std::fprintf(stderr, "event_phase=%zu\n", *result.event_phase);
            if (result.event_pair) std::fprintf(stderr, "event_pair=%zu,%zu\n", result.event_pair->first, result.event_pair->second);
            if (result.last_result) std::fprintf(stderr, "solver=%s raw_residual=%.17g classified_residual=%.17g\n",
                result.last_result->solver_status.c_str(), result.last_result->residual_norm(), result.last_result->convergence.residual);
            throw EquilibriumFailure(continuation_status(result.stop), observable);
        }
        context.central = result.states.back();
    }
    require_central(context, context.central, observable);
    context.compiled = std::make_unique<eq::CompiledSystem>(eq::compile_problem(*context.model, target));
    context.feed = feed;
    context.order = -1;
}

eq::SolvedStateActionBatch actions(EquilibriumContext& context, const std::vector<int>& coordinates) {
    eq::SolvedStateActionRequest request;
    for (std::size_t phase = 0; phase < 2; ++phase)
        for (std::size_t component = 0; component < context.map.offset.size(); ++component)
            request.observables.push_back({eq::SolvedStateObservableKind::PhaseComponentAmount, phase, component});
    for (const int coordinate : coordinates) {
        std::vector<double> direction(context.map.anchor.size(), 0.0);
        if (coordinate >= 0) direction.at(static_cast<std::size_t>(coordinate)) = 1.0;
        request.directions.push_back({0.0, 0.0, context.map.amounts(direction, true), {}});
    }
    auto batch = eq::solved_state_actions(*context.compiled, context.central, request);
    for (std::size_t row = 0; row < batch.results.size(); ++row) {
        const auto& item = batch.results[row];
        auto status = item.status;
        if (status == ActionStatus::Available && (!item.value || !item.action
            || !std::isfinite(*item.value) || !std::isfinite(*item.action))) status = ActionStatus::Nonfinite;
        if (status != ActionStatus::Available)
            throw EquilibriumFailure(status, eq::solved_state_observable_id(item.observable));
        context.values[row] = *item.value;
    }
    return batch;
}

double equilibrium_amount(arglist* arguments) {
    auto& prepared = prepared_equilibrium();
    std::lock_guard<std::mutex> lock(prepared.mutex);
    std::string observable = "request";
    if (arguments != nullptr) arguments->Errmsg = nullptr;
    invalidate(arguments);
    try {
        if (arguments == nullptr || arguments->nr < 1 || arguments->nr > 64 || arguments->n != arguments->nr + 3
            || arguments->at == nullptr || arguments->ra == nullptr || arguments->sa == nullptr)
            reject("equilibrium callback requires real controls and three strings");
        const int n = arguments->nr;
        for (int i = 0; i < arguments->n; ++i)
            if (arguments->at[i] != (i < n ? i : n - 1 - i)) reject("invalid equilibrium argument map");
        auto selector = [&](int index, std::size_t maximum) {
            const std::string s(string_argument(arguments, index));
            std::size_t v = 0;
            const auto parsed = std::from_chars(s.data(), s.data() + s.size(), v);
            if (s.empty() || parsed.ec != std::errc{} || parsed.ptr != s.data() + s.size() || v >= maximum)
                reject("invalid observable selector");
            return v;
        };
        const auto p = selector(n + 1, 2), c = selector(n + 2, 64);
        observable = eq::solved_state_observable_id({eq::SolvedStateObservableKind::PhaseComponentAmount, p, c});
        std::vector<double> feed(n);
        for (std::size_t i = 0; i < feed.size(); ++i) {
            feed[i] = arguments->ra[i];
            if (!std::isfinite(feed[i])) throw EquilibriumFailure(ActionStatus::Nonfinite, observable);
            if (feed[i] < 0.0) throw EquilibriumFailure(ActionStatus::OutsideInterval, observable);
        }
        const char* record = string_argument(arguments, n);
        if (!prepared.context) prepared.context = prepare_equilibrium(record, observable);
        auto& context = *prepared.context;
        if (context.record != record) reject("immutable equilibrium record changed");
        if (feed.size() != context.map.anchor.size() || c >= context.map.offset.size()) reject("control/observable dimension mismatch");
        if (context.order < 0 || context.feed != feed) select_feed(context, feed, observable);
        const int order = arguments->hes != nullptr ? 2 : (arguments->derivs != nullptr ? 1 : 0);
        if (order == 0 && context.order < 0) actions(context, {-1});
        if (order >= 1 && context.order < 1) {
            for (int i = 0; i < n; ++i) {
                const auto batch = actions(context, {i});
                for (std::size_t row = 0; row < context.values.size(); ++row) context.gradients[row][i] = *batch.results[row].action;
            }
        }
        if (order >= 2 && context.order < 2) {
            std::size_t packed = 0;
            for (int j = 0; j < n; ++j) for (int i = 0; i <= j; ++i) {
                const auto batch = actions(context, {i, j});
                for (std::size_t row = 0; row < context.values.size(); ++row) context.hessians[row][packed] = *batch.results[row].action;
                ++packed;
            }
        }
        context.order = std::max(context.order, order);
        const auto row = p * context.map.offset.size() + c;
        if (arguments->derivs != nullptr) std::copy(context.gradients[row].begin(), context.gradients[row].end(), arguments->derivs);
        if (arguments->hes != nullptr) std::copy(context.hessians[row].begin(), context.hessians[row].end(), arguments->hes);
        return context.values[row];  // Absolute phase-component mol on the caller's declared batch basis.
    } catch (const EquilibriumFailure& error) {
        callback_error = error.what();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "epcsaft_equilibrium diagnostic: %s\n", error.what());
        callback_error = EquilibriumFailure(ActionStatus::DefiningEquationsNotSatisfied, observable).what();
    } catch (...) {
        callback_error = EquilibriumFailure(ActionStatus::Nonfinite, observable).what();
    }
    // A solver only sees an evaluation error; name the typed cause and the trial controls.
    std::fprintf(stderr, "epcsaft_equilibrium failure: %s controls=", callback_error.c_str());
    for (int i = 0; arguments != nullptr && arguments->ra != nullptr && i < arguments->nr; ++i)
        std::fprintf(stderr, "%s%.17g", i != 0 ? "," : "", arguments->ra[i]);
    std::fputc('\n', stderr);
    if (prepared.context) prepared.context->order = -1;
    invalidate(arguments);
    if (arguments != nullptr) {
        if (arguments->hes != nullptr) callback_error.insert(0, "\"");
        else if (arguments->derivs != nullptr) callback_error.insert(0, "'");
        arguments->Errmsg = callback_error.data();
    }
    return invalid_value();
}

}  // namespace

extern "C" void funcadd_ASL(AmplExports* exports) {
    if (exports == nullptr || exports->Addfunc == nullptr) return;
    exports->Addfunc(
        "epcsaft_neutral_ar", residual_helmholtz,
        FUNCADD_REAL_VALUED | FUNCADD_STRING_ARGS, 6, nullptr, exports
    );
    exports->Addfunc(
        "epcsaft_neutral_pressure", pressure,
        FUNCADD_REAL_VALUED | FUNCADD_STRING_ARGS, 6, nullptr, exports
    );
    exports->Addfunc(
        "epcsaft_neutral_recordings", recordings,
        FUNCADD_REAL_VALUED, 0, nullptr, exports
    );
    exports->Addfunc(
        "epcsaft_equilibrium_phase_component_amount", equilibrium_amount,
        FUNCADD_REAL_VALUED | FUNCADD_STRING_ARGS, -5, nullptr, exports
    );
}
