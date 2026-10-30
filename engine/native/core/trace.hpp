#pragma once

#include "core/convergence.hpp"
#include "core/declared.hpp"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <initializer_list>
#include <limits>
#include <string>
#include <string_view>
#include <optional>
#include <utility>
#include <type_traits>
#include <vector>

namespace epcsaft::core {

struct Event final {
    const char* name;
    std::uint32_t span;
    std::uint32_t iteration;
    std::vector<double> values;
    std::string reason;
};

struct Span final {
    const char* name;
    std::optional<std::uint32_t> parent;
    std::optional<Convergence> convergence;
    double wall_us = 0.0;
    bool completed = false;
};

struct IpoptLog final {
    std::string destination, error;
    bool attempted = false;
    std::shared_ptr<void> application;
};

struct Trace final {
    std::vector<Span> spans;
    std::vector<Event> events;
    std::optional<std::uint32_t> current;
    bool diagnostics = false;
    bool incomplete = false;
    std::shared_ptr<IpoptLog> ipopt_log;
};

EPC_API Trace*& current_trace() noexcept;

inline void retain_first(std::optional<Convergence>& first, const std::optional<Convergence>& next, bool earlier = false) {
    if (next && (!first || earlier)) first = next;
}
template<class Factory, std::enable_if_t<std::is_invocable_v<Factory>, int> = 0>
void retain_first(std::optional<Convergence>& first, Factory&& build) {
    if (!first) first = build();
}
inline void retain_failure(std::optional<Convergence>& first, const std::exception& error, const char* operation) {
    if (first) return;
    Convergence record;
    record.residual = std::numeric_limits<double>::quiet_NaN();
    if (const auto* c = dynamic_cast<const ConvergenceFailure*>(&error)) record = c->first_numerical_failure.value_or(c->record);
    else if (const auto* d = dynamic_cast<const DeclarationError*>(&error); d && d->record) record = *d->record;
    if (record.cause.empty()) record.cause = std::string(operation) + ": " + error.what();
    first = std::move(record);
}
inline void retain_failure(std::optional<Convergence>& first, std::string_view cause, const char* operation) {
    if (!first) retain_failure(first, std::runtime_error(std::string(cause)), operation);
}
inline void retain_rank(std::vector<NumericalRank>& records, const NumericalRank& rank) {
    for (auto& r : records) if (r.basis == rank.basis) { r = rank; return; }
    records.push_back(rank);
}

template <typename Observation>
void observe(Trace* trace, Observation&& observation) {
    if (trace) observation();
}

class Capture final {
public:
    explicit Capture(std::optional<bool> diagnostics = std::nullopt, const std::string& ipopt_log = {}) : previous_(current_trace()) {
        if (!ipopt_log.empty()) { trace_.ipopt_log = std::make_shared<IpoptLog>(); trace_.ipopt_log->destination = ipopt_log; }
        else if (previous_) trace_.ipopt_log = previous_->ipopt_log;
        trace_.diagnostics = diagnostics ? *diagnostics : previous_ && previous_->diagnostics;
        parent_ = trace_.diagnostics && previous_ ? previous_->current : std::nullopt;
        current_trace() = &trace_;
    }
    ~Capture() noexcept(false) { close_log(); current_trace() = previous_; forward(); }
    Capture(const Capture&) = delete;
    Capture& operator=(const Capture&) = delete;

    const Trace& trace() const noexcept { return trace_; }
    Trace take() {
        close_log();
        forward();
        Trace result = std::move(trace_);
        trace_ = {};
        trace_.diagnostics = result.diagnostics;
        return result;
    }

private:
    Trace trace_;
    Trace* previous_;
    std::optional<std::uint32_t> parent_;
    bool forwarded_ = false;
    void close_log() noexcept {
        if (trace_.ipopt_log && (!previous_ || previous_->ipopt_log != trace_.ipopt_log)) trace_.ipopt_log->application.reset();
    }
    void forward() {
        if (forwarded_ || !previous_) return;
        forwarded_ = true;
        if (previous_->ipopt_log == trace_.ipopt_log) previous_->incomplete |= trace_.incomplete;
        if (!trace_.diagnostics || !previous_->diagnostics) return;
        const auto spans = previous_->spans.size();
        for (auto span : trace_.spans) {
            span.parent = span.parent ? std::optional<std::uint32_t>(*span.parent + spans) : parent_;
            previous_->spans.push_back(std::move(span));
        }
        for (auto event : trace_.events) {
            event.span += spans;
            previous_->events.push_back(std::move(event));
        }
    }
};

class Scope final {
public:
    explicit Scope(const char* name) : trace_(current_trace()) {
        if (trace_ == nullptr) return;
        const auto now = std::chrono::steady_clock::now();
        parent_ = trace_->current;
        observe([&] { trace_->spans.push_back({name, parent_, std::nullopt, 0.0}); });
        index_ = static_cast<std::uint32_t>(trace_->spans.size() - 1);
        trace_->current = index_;
        start_ = now;
    }

    ~Scope() {
        if (trace_ == nullptr) return;
        const auto now = std::chrono::steady_clock::now();
        trace_->spans[*index_].wall_us =
            std::chrono::duration<double, std::micro>(now - start_).count();
        trace_->spans[*index_].completed = std::uncaught_exceptions() == exceptions_;
        trace_->current = parent_;
    }

    bool on() const noexcept { return trace_ != nullptr; }
    bool diagnostics_on() const noexcept { return on() && trace_->diagnostics; }
    static std::optional<Scope> detailed(const char* name) {
        const auto* trace = current_trace();
        if (!trace || !trace->diagnostics) return std::nullopt;
        return std::optional<Scope>(std::in_place, name);
    }
    bool nested(const char* name) const noexcept {
        for (auto parent = parent_; parent; parent = trace_->spans[*parent].parent)
            if (std::strcmp(trace_->spans[*parent].name, name) == 0) return true;
        return false;
    }
    template <typename Observation>
    void observe(Observation&& observation) {
        core::observe(trace_, std::forward<Observation>(observation));
    }

    void probe(const char* name, double value, std::uint32_t iteration = 0) {
        probe(name, {value}, iteration);
    }
    void probe(const char* name, std::initializer_list<double> values,
               std::uint32_t iteration = 0, const char* reason = "") {
        observe([&] { trace_->events.push_back({name, *index_, iteration, values, reason}); });
    }

    void probe(const char* name, const double* values, std::size_t count,
               std::uint32_t iteration = 0) {
        observe([&] { trace_->events.push_back({
            name, *index_, iteration, std::vector<double>(values, values + count), {}
        }); });
    }

    void probe(const char* name, const std::vector<double>& values,
               std::uint32_t iteration = 0) {
        observe([&] { trace_->events.push_back({name, *index_, iteration, values, {}}); });
    }

    void finish(const Convergence& record) {
        observe([&] { trace_->spans[*index_].convergence = record; });
    }

private:
    Trace* trace_;
    std::optional<std::uint32_t> index_;
    std::optional<std::uint32_t> parent_;
    std::chrono::steady_clock::time_point start_;
    int exceptions_ = std::uncaught_exceptions();
};

}  // namespace epcsaft::core
