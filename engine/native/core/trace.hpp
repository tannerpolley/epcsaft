#pragma once

#include "core/convergence.hpp"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace epcsaft::core {

struct Event final {
    const char* name;
    std::uint32_t span;
    std::uint32_t iteration;
    std::vector<double> values;
};

struct Span final {
    const char* name;
    std::optional<std::uint32_t> parent;
    std::optional<Convergence> convergence;
    double wall_us = 0.0;
};

struct Trace final {
    std::vector<Span> spans;
    std::vector<Event> events;
    std::optional<std::uint32_t> current;
};

inline thread_local Trace* active_trace = nullptr;

class Capture final {
public:
    Capture() : previous_(active_trace) { active_trace = &trace_; }
    ~Capture() { active_trace = previous_; }
    Capture(const Capture&) = delete;
    Capture& operator=(const Capture&) = delete;

    const Trace& trace() const noexcept { return trace_; }
    Trace take() {
        Trace result = std::move(trace_);
        trace_ = {};
        return result;
    }

private:
    Trace trace_;
    Trace* previous_;
};

class Scope final {
public:
    explicit Scope(const char* name) noexcept : trace_(active_trace) {
        if (trace_ == nullptr) return;
        const auto now = std::chrono::steady_clock::now();
        parent_ = trace_->current;
        trace_->spans.push_back({name, parent_, std::nullopt, 0.0});
        index_ = static_cast<std::uint32_t>(trace_->spans.size() - 1);
        trace_->current = index_;
        start_ = now;
    }

    ~Scope() {
        if (trace_ == nullptr) return;
        const auto now = std::chrono::steady_clock::now();
        trace_->spans[*index_].wall_us =
            std::chrono::duration<double, std::micro>(now - start_).count();
        trace_->current = parent_;
    }

    bool on() const noexcept { return trace_ != nullptr; }

    void probe(const char* name, double value, std::uint32_t iteration = 0) {
        if (trace_ == nullptr) return;
        trace_->events.push_back({name, *index_, iteration, {value}});
    }

    void probe(const char* name, const double* values, std::size_t count,
               std::uint32_t iteration = 0) {
        if (trace_ == nullptr) return;
        trace_->events.push_back({
            name, *index_, iteration, std::vector<double>(values, values + count)
        });
    }

    void probe(const char* name, const std::vector<double>& values,
               std::uint32_t iteration = 0) {
        if (trace_ == nullptr) return;
        trace_->events.push_back({name, *index_, iteration, values});
    }

    void finish(const Convergence& record) {
        if (trace_ == nullptr) return;
        trace_->spans[*index_].convergence = record;
    }

private:
    Trace* trace_;
    std::optional<std::uint32_t> index_;
    std::optional<std::uint32_t> parent_;
    std::chrono::steady_clock::time_point start_;
};

}  // namespace epcsaft::core
