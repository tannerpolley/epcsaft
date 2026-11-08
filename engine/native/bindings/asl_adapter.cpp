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

using epcsaft::DerivativeOrder;
using epcsaft::HelmholtzModel;
using epcsaft::NeutralProperty;
using epcsaft::NeutralPropertyDerivatives;
using epcsaft::detail::family_values;

struct PreparedModel final {
    std::mutex mutex;
    std::tuple<std::string, std::string, std::string, std::string> inputs;
    std::unique_ptr<HelmholtzModel> model;
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

double real_argument(arglist* arguments, int position) {
    const int index = arguments->at[position];
    if (index < 0 || index >= arguments->nr) reject("invalid real argument map");
    return arguments->ra[index];
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
int argument_count(arglist* arguments, int strings, int maximum, const char* message) {
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

HelmholtzModel& model_for(
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
    prepared.model = std::make_unique<HelmholtzModel>(
        epcsaft::build_parameters(document)
    );
    prepared.inputs = inputs;
    return *prepared.model;
}

double evaluate(
    arglist* arguments,
    NeutralProperty property
) {
    const char* message = "epcsaft neutral function requires two real and four string arguments";
    if (argument_count(arguments, 4, 2, message) != 2) reject(message);
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
    template<class Read> auto optional(Read read) {  // a vector of length 0 or 1
        const auto v = vector(read, 1);
        return v.empty() ? std::nullopt : std::optional(v.front());
    }
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
    // The record's model at theta0, which owns anchor_result, and the model at another theta.
    std::unique_ptr<HelmholtzModel> model, shifted;
    eq::EquilibriumProblem anchor;
    eq::ContinuationOptions options;
    eq::AslFeedMap map;
    eq::EquilibriumResult anchor_result;
    std::unique_ptr<eq::CompiledSystem> compiled;
    eq::EquilibriumResult central;
    std::vector<double> feed;  // the q of `compiled` and `central`; a null `compiled` is no valid solved state
    struct Rows final {
        std::vector<double> values;
        std::vector<std::vector<double>> gradients, hessians;
        int order = -1;  // highest order filled at `feed`
    } rows[2];  // 0: the 2C phase-component amounts, 1: the system enthalpy
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

void require_central(const EquilibriumContext& context, const eq::EquilibriumResult& result, const std::string& observable) {
    const auto failure = [&](ActionStatus status) { return EquilibriumFailure(status, observable, result.message + conditioning_warnings(result.rank_records)
        + (result.first_numerical_failure ? "; first refusal: " + epcsaft::core::describe(*result.first_numerical_failure) : "")); };
    const auto count = context.map.offset.size(), identity = context.map.identity_component;
    if (result.topology_event == eq::TopologyEvent::PhasesCoalesced)
        throw failure(ActionStatus::NonuniqueCoalescedPhases);
    if (result.topology_event == eq::TopologyEvent::PhaseDisappeared)
        throw failure(ActionStatus::PhaseBoundaryNotDifferentiable);
    if (!result.has_state()) throw failure(ActionStatus::CentralSolveUnavailable);
    if (!result.requested_tolerance_met.value_or(false))
        throw failure(ActionStatus::DefiningEquationsNotSatisfied);
    if (result.topology_event != eq::TopologyEvent::None || result.phase_count != 2
        || result.mole_fractions.size() != 2 * count || result.amounts.size() != 2
        || !(result.mole_fractions[identity] - result.mole_fractions[count + identity] > context.map.minimum_gap))
        throw failure(ActionStatus::PhaseBoundaryNotDifferentiable);
}

epcsaft::ParameterDocument read_document(RecordReader& reader) {
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
    reader.field("polar_choice");
    document.polar_choice = reader.next();
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
    return document;
}

std::vector<epcsaft::ActiveParameter> read_active(RecordReader& reader) {
    reader.field("active_parameters");
    return reader.vector([&] {
        const auto family = epcsaft::ParameterFamily(reader.integer(epcsaft::detail::PARAMETER_FAMILY_COUNT - 1));
        return epcsaft::ActiveParameter{family, reader.vector([&] { return reader.integer(SIZE_MAX); })};
    });
}

// The writers' put_reaction and put_source_reference records.
eq::ReactionSpec read_reaction(RecordReader& reader) {
    eq::ReactionSpec r;
    r.name = reader.next();
    r.stoichiometry = reader.vector([&] { return reader.number(); });
    const auto target = reader.next();
    if (target == "K") r.equilibrium_constant = reader.number();
    else if (target == "affinity") r.affinity_j_per_mol = reader.number();
    else if (target == "correlation") {
        auto& c = r.correlation.emplace();
        for (double* v : eq::law_fields(c)) *v = reader.number();
        c.standard_state_id = reader.next();
    } else reject("unsupported reaction target");
    if (reader.integer(1) != 0) r.reference = {reader.optional([&] { return reader.number(); }), epcsaft::ReferenceSourceBasis(reader.integer(1))};
    return r;
}

std::optional<epcsaft::NeutralReferenceDeclaration> read_source_reference(RecordReader& reader) {
    reader.field("source_reference");
    if (reader.integer(1) == 0) return std::nullopt;  // braced initialization reads left to right
    return epcsaft::NeutralReferenceDeclaration{reader.optional([&] { return reader.integer(); }), reader.vector([&] { return reader.number(); }),
        reader.number(), reader.vector([&] { return epcsaft::SpeciesReference(reader.integer(2)); })};
}

// The block put_thermochemistry writes.
std::optional<epcsaft::ThermochemistryRecord> read_thermochemistry(RecordReader& reader) {
    if (reader.integer(1) == 0) return std::nullopt;
    epcsaft::ThermochemistryRecord thermochemistry;
    reader.field("standard_pressure");
    thermochemistry.standard_pressure = reader.number();
    reader.field("component_intervals");
    thermochemistry.component_intervals = reader.vector([&] {
        return reader.vector([&] {
            epcsaft::IdealInterval item{reader.number(), reader.number(), reader.integer(1) != 0, reader.integer(1) != 0, {}};
            item.correlation.source_pressure = reader.number();
            const auto tag = reader.next();
            const auto form = std::find(eq::ideal_tags.begin(), eq::ideal_tags.end(), tag) - eq::ideal_tags.begin();
            if (form == 5) reject("unknown ASL thermochemistry form " + tag);
            item.correlation.form = eq::ideal_forms[form];
            std::visit([&](auto& c) {
                const auto count = reader.integer();
                if constexpr (std::is_same_v<decltype(c), epcsaft::IdealPolynomial&>) c.coefficients.resize(count);
                if (count != c.coefficients.size()) reject("invalid ASL interval");
                for (double& value : c.coefficients) value = reader.number();
                if (reader.integer() != eq::ideal_anchors(c).size()) reject("invalid ASL interval");
                for (double* value : eq::ideal_anchors(c)) *value = reader.number();
            }, item.correlation.form);
            return item;
        });
    });
    return thermochemistry;
}

// The model at theta, by the slot assignment of CoordinateSet::apply, in model.active_parameters() order.
auto model_at(const HelmholtzModel& model, const std::vector<double>& theta) {
    auto parameters = model.parameters();
    const auto& active = model.active_parameters();
    for (std::size_t k = 0; k < theta.size(); ++k)
        for (auto i : active[k].indices) family_values(parameters, active[k].family).first[i] = theta[k];
    return std::make_unique<HelmholtzModel>(std::move(parameters), active, model.thermochemistry(), model.solver_options());
}

std::unique_ptr<EquilibriumContext> prepare_equilibrium(const char* raw, const std::string& observable) {
    RecordReader reader(raw);
    reader.field("epcsaft.asl.equilibrium.v7");
    const auto document = read_document(reader);
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
        const auto kind = reader.next();
        if (kind != "liquid" && kind != "vapor") reject("equilibrium callback requires liquid or vapor phases");
        phase.kind = kind == "liquid" ? eq::PhaseKind::Liquid : eq::PhaseKind::Vapor;
        phase.support = reader.vector([&] { return reader.integer(); });
        phase.support_declared = reader.integer(1) != 0;
        if (reader.next() != "free") reject("equilibrium callback requires free phase amounts");
        phase.amount.fraction = reader.number();
        phase.composition_guess = reader.vector([&] { return reader.number(); });
        phase.packing_guess = reader.optional([&] { return reader.number(); });
        return phase;
    });
    reader.field("reactions");
    problem.reactions = reader.vector([&] { return read_reaction(reader); });
    reader.field("standard_gibbs_j_per_mol");
    problem.standard_gibbs_j_per_mol = reader.vector([&] { return reader.number(); });
    problem.source_reference = read_source_reference(reader);
    reader.field("topology_declared");
    problem.topology_declared = reader.integer(1) != 0;
    reader.field("solver_options");
    for (auto& value : problem.solver_options.values) value = reader.number();
    const auto active = read_active(reader);
    const auto thermochemistry = read_thermochemistry(reader);
    reader.field("continuation_options");
    context->options.initial_step = reader.number();
    context->options.minimum_step = reader.number();
    context->options.max_steps = reader.integer(SIZE_MAX);  // continue_equilibrium validates both counts
    context->options.max_retries_per_step = reader.integer(SIZE_MAX);
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
    context->model = std::make_unique<HelmholtzModel>(epcsaft::build_parameters(document), active, thermochemistry, problem.solver_options);
    context->map.validate(problem, context->model->parameters());
    const std::size_t counts[2] = {2 * document.component_ids.size(), 1}, n = 2 + context->map.anchor.size() + active.size();
    for (int kind = 0; kind < 2; ++kind) {
        auto& rows = context->rows[kind];
        rows.values.resize(counts[kind]);
        rows.gradients.assign(counts[kind], std::vector<double>(n));
        rows.hessians.assign(counts[kind], std::vector<double>(n * (n + 1) / 2));
    }
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
// anchor (T, P, feed) precedes one leg moving T, P and the feed together, since a parameter path moves nothing else.
void select_state(EquilibriumContext& context, const std::vector<double>& q, const std::string& observable) {
    const auto controls = context.map.anchor.size();
    const std::vector<double> feed(q.begin() + 2, q.begin() + 2 + controls), theta(q.begin() + 2 + controls, q.end());
    context.compiled.reset();  // it points at its model
    if (theta == context.model->active_parameter_values()) context.shifted.reset();
    if (theta != (context.shifted ? *context.shifted : *context.model).active_parameter_values()) try {
        context.shifted = model_at(*context.model, theta);  // rebuilt against the current model's own values
    } catch (const std::invalid_argument& error) { throw EquilibriumFailure(ActionStatus::OutsideInterval, observable, error.what()); }
    const auto& model = context.shifted ? *context.shifted : *context.model;
    auto target = context.anchor;
    target.temperature = eq::FreeValue::fixed(q[0]), target.pressure = eq::FreeValue::fixed(q[1]);
    target.feed.values = context.map.amounts(feed);
    for (const double amount : target.feed.values)
        if (!std::isfinite(amount) || amount < 0.0) throw EquilibriumFailure(ActionStatus::OutsideInterval, observable);
    const bool moved = feed != context.map.anchor || q[0] != context.anchor.temperature.value || q[1] != context.anchor.pressure.value;
    if (!moved && !context.shifted) context.central = context.anchor_result;
    else try {
        context.central = eq::solve_equilibrium(model, target);
        require_central(context, context.central, observable);
    } catch (const std::exception& original) {
        const auto leg = [&](const auto& path, const auto& to, const auto& targets, const char* name) {
            context.central = continuation_leg(context, path, to, targets, original.what(), name, observable);
        };
        context.central = context.anchor_result;
        if (context.shifted) leg(*context.model, context.anchor, theta, "parameter continuation");
        if (moved) leg(model, target, std::vector<double>{}, "continuation");
    }
    require_central(context, context.central, observable);
    context.compiled = std::make_unique<eq::CompiledSystem>(eq::compile_problem(model, target));
    context.feed = q;
    for (auto& rows : context.rows) rows.order = -1;
}

// One native batch of one kind's rows; an unavailable row fails only its own kind's callback.
// ponytail: a caller reading both kinds at a new q pays both batches (185 ms Hessian against 116 ms combined); the upgrade is
// one combined batch with per-row failure orders, as PropertyContext::Branch::failures keeps, which could also fill the
// amount rows from an enthalpy call when they are all available.
eq::SolvedStateActionBatch actions(EquilibriumContext& context, const std::vector<int>& coordinates, int kind) {
    eq::SolvedStateActionRequest request;
    if (kind == 1) request.observables.push_back({eq::SolvedStateObservableKind::SystemEnthalpy});
    for (std::size_t phase = 0; kind == 0 && phase < 2; ++phase)
        for (std::size_t component = 0; component < context.map.offset.size(); ++component)
            request.observables.push_back({eq::SolvedStateObservableKind::PhaseComponentAmount, phase, component});
    const auto controls = context.map.anchor.size();
    for (const int coordinate : coordinates) {  // a q column (T, P, feed controls, theta); -1 is no direction
        std::vector<double> e(2 + controls + context.model->active_parameters().size(), 0.0);
        if (coordinate >= 0) e.at(coordinate) = 1.0;
        const auto theta = e.begin() + 2 + controls;
        request.directions.push_back({e[0], e[1], context.map.amounts({e.begin() + 2, theta}, true), {theta, e.end()}});
    }
    auto batch = eq::solved_state_actions(*context.compiled, context.central, request);
    for (std::size_t row = 0; row < batch.results.size(); ++row) {
        const auto& item = batch.results[row];
        auto status = item.status;
        if (status == ActionStatus::Available && (!item.value || !item.action
            || !std::isfinite(*item.value) || !std::isfinite(*item.action))) status = ActionStatus::Nonfinite;
        if (status != ActionStatus::Available) {
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
    return invalid_value();
}

// kind 0: (q, record, phase, component), mol; kind 1: (q, record), the system enthalpy in J. Both kinds share one
// solved state per record and q.
double equilibrium_observable(arglist* arguments, int kind) {
    static Prepared<EquilibriumContext> prepared;
    std::lock_guard<std::mutex> lock(prepared.mutex);
    std::string observable = "request";
    if (arguments != nullptr) arguments->Errmsg = nullptr;
    invalidate(arguments);
    try {
        const int n = argument_count(arguments, kind == 0 ? 3 : 1, 128, "equilibrium callback requires real q and the record");
        const auto p = kind == 0 ? selector(arguments, n + 1, 2) : 0, c = kind == 0 ? selector(arguments, n + 2, 64) : 0;
        observable = eq::solved_state_observable_id(kind == 0 ? eq::SolvedStateObservable{eq::SolvedStateObservableKind::PhaseComponentAmount, p, c}
            : eq::SolvedStateObservable{eq::SolvedStateObservableKind::SystemEnthalpy});
        const char* record = string_argument(arguments, n);
        auto& context = prepared.find(record, "equilibrium", [&] { return prepare_equilibrium(record, observable); });
        const auto controls = 2 + context.map.anchor.size();  // T, P and the feed controls
        if (std::size_t(n) != controls + context.model->active_parameters().size() || c >= context.map.offset.size())
            reject("control/observable dimension mismatch");
        const std::vector<double> q(arguments->ra, arguments->ra + n);  // a parameter may be negative (k_ij)
        for (std::size_t i = 0; i < q.size(); ++i) {
            if (!std::isfinite(q[i])) throw EquilibriumFailure(ActionStatus::Nonfinite, observable);
            if (i < controls && (i < 2 ? q[i] <= 0.0 : q[i] < 0.0)) throw EquilibriumFailure(ActionStatus::OutsideInterval, observable);
        }
        if (!context.compiled || context.feed != q) select_state(context, q, observable);
        auto& rows = context.rows[kind];
        const int order = arguments->hes != nullptr ? 2 : (arguments->derivs != nullptr ? 1 : 0);
        if (order == 0 && rows.order < 0) actions(context, {-1}, kind);
        if (order >= 1 && rows.order < 1) {
            for (int i = 0; i < n; ++i) {
                const auto batch = actions(context, {i}, kind);
                for (std::size_t row = 0; row < rows.values.size(); ++row) rows.gradients[row][i] = *batch.results[row].action;
            }
        }
        if (order >= 2 && rows.order < 2) {
            std::size_t packed = 0;
            for (int j = 0; j < n; ++j) for (int i = 0; i <= j; ++i) {
                const auto batch = actions(context, {i, j}, kind);
                for (std::size_t row = 0; row < rows.values.size(); ++row) rows.hessians[row][packed] = *batch.results[row].action;
                ++packed;
            }
        }
        rows.order = std::max(rows.order, order);
        const auto row = p * context.map.offset.size() + c;
        if (arguments->derivs != nullptr) std::copy(rows.gradients[row].begin(), rows.gradients[row].end(), arguments->derivs);
        if (arguments->hes != nullptr) std::copy(rows.hessians[row].begin(), rows.hessians[row].end(), arguments->hes);
        return rows.values[row];  // On the caller's declared absolute batch basis.
    } catch (const EquilibriumFailure& error) {
        callback_error = error.what();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "epcsaft_equilibrium diagnostic: %s\n", error.what());
        callback_error = EquilibriumFailure(ActionStatus::DefiningEquationsNotSatisfied, observable, error.what()).what();
    } catch (...) {
        callback_error = EquilibriumFailure(ActionStatus::Nonfinite, observable).what();
    }
    for (const auto& item : prepared.contexts) item->compiled.reset();  // every record re-solves at its next call
    return report_failure(arguments, "epcsaft_equilibrium");
}

double equilibrium_amount(arglist* arguments) { return equilibrium_observable(arguments, 0); }
double equilibrium_enthalpy(arglist* arguments) { return equilibrium_observable(arguments, 1); }

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
    } branches[2], solid_rows;  // liquid (highest stable root), vapor (lowest); solid rows need no density root
    std::vector<eq::SolidSpec> solids;
    std::optional<epcsaft::NeutralReferenceDeclaration> reference;
};

std::unique_ptr<PropertyContext> prepare_property(const char* raw) {
    RecordReader reader(raw);
    reader.field("epcsaft.asl.property.v5");
    const bool require_phase_identity = reader.integer(1) != 0;
    const auto document = read_document(reader);
    epcsaft::policy::SolverOptions solver_options;
    reader.field("solver_options");
    for (auto& value : solver_options.values) value = reader.number();
    const auto active = read_active(reader);
    const auto thermochemistry = read_thermochemistry(reader);
    auto context = std::make_unique<PropertyContext>();
    reader.field("solids");
    const auto count = reader.integer();
    if (count > 0) context->reference = read_source_reference(reader);
    for (std::size_t s = 0; s < count; ++s) {
        const auto r = read_reaction(reader);
        if (!r.correlation) reject("ASL solid rows require a ln K(T) law");
        context->solids.push_back({r.name, r.stoichiometry, r.equilibrium_constant, r.correlation, r.reference, reader.number(), {}});
    }
    reader.finish();
    context->require_phase_identity = require_phase_identity;
    context->model = std::make_unique<HelmholtzModel>(
        epcsaft::build_parameters(document), active, thermochemistry, solver_options);
    eq::admit_solid_declarations(*context->model, context->reference, context->solids);
    context->record = raw;
    return context;
}

void evaluate_property(PropertyContext& context, std::size_t phase, const std::vector<double>& q, int order) {
    auto& branch = context.branches[phase];
    if (branch.q == q && branch.order >= order) return;
    using Status = epcsaft::PropertyActionStatus;
    const std::size_t c = context.model->parameters().component_count(), rows = 2 + c, n = q.size();
    const auto active = context.model->active_parameters();
    const std::vector<double> theta(q.begin() + rows, q.end());
    if (theta != context.model->active_parameter_values()) context.model = model_at(*context.model, theta);
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
    branch.q = q;
    branch.order = -1;
    branch.values.assign(rows, 0.0);
    branch.failures.assign(rows, {3, std::string()});
    branch.gradients.assign(rows, std::vector<double>(n));
    branch.hessians.assign(rows, std::vector<double>(n * (n + 1) / 2));
    const auto run = [&](const std::vector<epcsaft::PressureActionDirection>& directions) {
        const auto batch = context.model->pressure_actions(q[0], q[1], x,
            phase == 0 ? epcsaft::DensityBranch::Liquid : epcsaft::DensityBranch::Vapor, observables, directions, anchor);
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
        const auto requested = phase == 0 ? epcsaft::DensityBranch::Liquid : epcsaft::DensityBranch::Vapor;
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
        const bool first = !anchor;
        if (first) anchor = branch.values[0] = *batch.molar_density;
        actions[0] = directions.size() == 1 ? batch.density_action.at(0) : batch.density_mixed_action.value();
        for (std::size_t row = 1; row < rows; ++row) {
            const auto& item = batch.results[row - 1];
            if (first && item.value) branch.values[row] = *item.value;  // an owned thermochemistry endpoint keeps its value
            if (first && !item.value) fail(row, 0, item.status);
            actions[row] = item.action.value_or(0.0);
            if (item.status != Status::Available || !item.action) fail(row, static_cast<int>(directions.size()), item.status);
        }
        return actions;
    };
    for (std::size_t k = 0; k < (order >= 1 ? n : 1); ++k) {
        const auto actions = run({basis[k]});
        for (std::size_t row = 0; row < rows; ++row) branch.gradients[row][k] = actions[row] / std::pow(sum, k >= 2 && k < rows);
    }
    std::size_t packed = 0;
    for (std::size_t j = 0; order >= 2 && j < n; ++j)
        for (std::size_t i = 0; i <= j; ++i, ++packed) {
            const auto actions = run({basis[i], basis[j]});
            for (std::size_t row = 0; row < rows; ++row)
                branch.hessians[row][packed] = actions[row] / std::pow(sum, (i >= 2 && i < rows) + (j >= 2 && j < rows))
                    - (i >= 2 && j < rows ? (sigma[j - 2] * branch.gradients[row][i] + sigma[i - 2] * branch.gradients[row][j]) / sum : 0.0);
        }
    branch.order = order;
}

// EqID solid_standard_gibbs on the record basis: l_s = ln K^rho_s(T, P0) + V_s (P - P0)/RT - dnu ln(P/(RT rho0)), rho0 =
// 1 mol/m^3, so sum_i nu_i (ln x_i + ln phi_i) = l_s at saturation. It depends on (T, P, theta) only.
void evaluate_solids(PropertyContext& context, const std::vector<double>& q, int order) {
    auto& rows = context.solid_rows;
    if (rows.q == q && rows.order >= order) return;
    const std::size_t S = context.solids.size(), n = q.size(), first = n - context.model->active_parameters().size();
    const std::vector<double> theta(q.begin() + first, q.end());
    if (theta != context.model->active_parameter_values()) context.model = model_at(*context.model, theta);
    std::vector<std::size_t> pairs(order >= 2 ? n - first : 0);
    std::iota(pairs.begin(), pairs.end(), 0);
    rows = {q, std::vector<double>(S), std::vector<std::vector<double>>(S, std::vector<double>(n)),
            std::vector<std::vector<double>>(S, std::vector<double>(n * (n + 1) / 2)), std::vector<std::pair<int, std::string>>(S, {3, ""}), order};
    const double T = q[0], P = q[1], RT = epcsaft::GAS_CONSTANT_J_PER_MOL_K * T;
    for (std::size_t s = 0; s < S; ++s) {
        const auto& solid = context.solids[s];
        const auto row = std::to_string(first + s) + ": solid " + solid.name + ": ";
        const auto fail = [&](int k, const char* cause, const std::string& detail) {
            if (k < rows.failures[s].first) rows.failures[s] = {k, "epcsaft_property:" + (cause + (":" + row)) + detail};
        };
        const double P0 = *solid.reference->reference_pressure_pa;
        eq::SolidLogK law;
        try {  // the reference is evaluated at its fixed source pressure, so it is independent of P
            law = eq::solid_log_k(*context.model, T, P0, context.reference, {solid}, std::max(order, 0), order >= 1 && first < n, pairs)[0];
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
            const auto v = C.derivative_value(a, 0, i, j);
            if (!v) fail(k, "unavailable_derivative", "reference column d^" + std::to_string(a) + "/dT; " + C.message);
            return v.value_or(0.0);
        };
        auto& g = rows.gradients[s];
        const auto h = [&](std::size_t i, std::size_t j) -> double& { return rows.hessians[s][j * (j + 1) / 2 + i]; };  // i <= j
        if (order >= 1) {
            g[0] = *law.slope + column(1, 1) - V * (P - P0) / (RT * T) + dnu / T;
            g[1] = V / RT - dnu / P;
            for (std::size_t k = first; k < n; ++k) g[k] = column(1, 0, k - first);
        }
        if (order < 2) continue;
        h(0, 0) = *law.curvature + column(2, 2) + 2.0 * V * (P - P0) / (RT * T * T) - dnu / (T * T);
        h(0, 1) = -V / (RT * T);
        h(1, 1) = dnu / (P * P);
        for (std::size_t k = first; k < n; ++k) {
            h(0, k) = column(2, 1, k - first);
            for (std::size_t l = first; l <= k; ++l) h(l, k) = column(2, 0, l - first, k - first);
        }
    }
}

double pressure_property(arglist* arguments) {
    static Prepared<PropertyContext> prepared;
    std::lock_guard<std::mutex> lock(prepared.mutex);
    if (arguments != nullptr) arguments->Errmsg = nullptr;
    invalidate(arguments);
    try {
        const int n = argument_count(arguments, 3, 130, "property callback requires real arguments and three strings");
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
        const int order = arguments->hes != nullptr ? 2 : (arguments->derivs != nullptr ? 1 : 0);
        const bool solid = row >= 2 + c;
        solid ? evaluate_solids(context, q, order) : evaluate_property(context, phase, q, order);
        const auto& branch = solid ? context.solid_rows : context.branches[phase];
        const auto r = solid ? row - 2 - c : row;
        if (branch.failures[r].first <= order) throw std::runtime_error(branch.failures[r].second);
        if (arguments->derivs != nullptr) std::copy(branch.gradients[r].begin(), branch.gradients[r].end(), arguments->derivs);
        if (arguments->hes != nullptr) std::copy(branch.hessians[r].begin(), branch.hessians[r].end(), arguments->hes);
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
    exports->Addfunc(
        "epcsaft_equilibrium_system_enthalpy", equilibrium_enthalpy,
        FUNCADD_REAL_VALUED | FUNCADD_STRING_ARGS, -3, nullptr, exports
    );
    exports->Addfunc(
        "epcsaft_pressure_property", pressure_property,
        FUNCADD_REAL_VALUED | FUNCADD_STRING_ARGS, -7, nullptr, exports
    );
}
