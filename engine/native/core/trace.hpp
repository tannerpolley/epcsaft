#pragma once

#include "core/convergence.hpp"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <initializer_list>
#include <string>
#include <optional>
#include <utility>
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

struct Trace final {
    std::vector<Span> spans;
    std::vector<Event> events;
    std::optional<std::uint32_t> current;
    bool diagnostics = false;
    bool incomplete = false;
};

Trace*& current_trace() noexcept;

template <typename Observation>
void observe(Trace* trace, Observation&& observation) noexcept {
    if (trace == nullptr) return;
    try { observation(); } catch (...) { trace->incomplete = true; }
}

class Capture final {
public:
    explicit Capture(std::optional<bool> diagnostics = std::nullopt) noexcept : previous_(current_trace()) {
        trace_.diagnostics = diagnostics ? *diagnostics : previous_ && previous_->diagnostics;
        parent_ = trace_.diagnostics && previous_ ? previous_->current : std::nullopt;
        current_trace() = &trace_;
    }
    ~Capture() { forward(); current_trace() = previous_; }
    Capture(const Capture&) = delete;
    Capture& operator=(const Capture&) = delete;

    const Trace& trace() const noexcept { return trace_; }
    Trace take() noexcept {
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
    void forward() noexcept {
        if (forwarded_ || !trace_.diagnostics || !previous_ || !previous_->diagnostics) return;
        forwarded_ = true;
        const auto spans = previous_->spans.size(), events = previous_->events.size();
        observe(previous_, [&] {
            previous_->incomplete = previous_->incomplete || trace_.incomplete;
            try {
                for (auto span : trace_.spans) {
                    span.parent = span.parent ? std::optional<std::uint32_t>(*span.parent + spans) : parent_;
                    previous_->spans.push_back(std::move(span));
                }
                for (auto event : trace_.events) {
                    event.span += spans;
                    previous_->events.push_back(std::move(event));
                }
            } catch (...) {
                previous_->spans.resize(spans); previous_->events.resize(events);
                trace_.incomplete = true;
                throw;
            }
        });
    }
};

class Scope final {
public:
    explicit Scope(const char* name) noexcept : trace_(current_trace()) {
        if (trace_ == nullptr) return;
        const auto now = std::chrono::steady_clock::now();
        parent_ = trace_->current;
        const auto size = trace_->spans.size();
        observe([&] { trace_->spans.push_back({name, parent_, std::nullopt, 0.0}); });
        if (trace_->spans.size() == size) { trace_ = nullptr; return; }
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
    static std::optional<Scope> detailed(const char* name) noexcept {
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
    void observe(Observation&& observation) noexcept {
        core::observe(trace_, std::forward<Observation>(observation));
    }

    void probe(const char* name, double value, std::uint32_t iteration = 0) noexcept {
        probe(name, {value}, iteration);
    }
    void probe(const char* name, std::initializer_list<double> values,
               std::uint32_t iteration = 0, const char* reason = "") noexcept {
        observe([&] { trace_->events.push_back({name, *index_, iteration, values, reason}); });
    }

    void probe(const char* name, const double* values, std::size_t count,
               std::uint32_t iteration = 0) noexcept {
        observe([&] { trace_->events.push_back({
            name, *index_, iteration, std::vector<double>(values, values + count), {}
        }); });
    }

    void probe(const char* name, const std::vector<double>& values,
               std::uint32_t iteration = 0) noexcept {
        observe([&] { trace_->events.push_back({name, *index_, iteration, values, {}}); });
    }

    void finish(const Convergence& record) noexcept {
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
