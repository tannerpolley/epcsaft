#pragma once
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iomanip>
#include <initializer_list>
#include <limits>
#include <type_traits>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace ep235 {
using Clock = std::chrono::steady_clock;
inline std::int64_t now() {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now().time_since_epoch()).count();
}
inline thread_local bool active = false, tracing = false;
inline thread_local unsigned stage = 0, trial = 0, trial_sequence = 0;
inline thread_local std::uint64_t sequence = 0;
inline thread_local char family = '-';
inline thread_local std::vector<std::string> events;
struct Key {
    const char* name;
    unsigned stage, first, second, third;
    char family;
    int outcome;
    bool operator==(const Key& b) const {
        return name == b.name && stage == b.stage && first == b.first && second == b.second
            && third == b.third && family == b.family && outcome == b.outcome;
    }
};
struct Hash {
    std::size_t operator()(const Key& k) const {
        std::size_t h = std::hash<const char*>{}(k.name);
        for (unsigned v : {k.stage, k.first, k.second, k.third, unsigned(k.family), unsigned(k.outcome)})
            h ^= v + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};
struct Cell {
    std::uint64_t count = 0;
    std::int64_t exclusive = 0, inclusive = 0;
};
inline thread_local std::unordered_map<Key, Cell, Hash> cells;

inline void quoted(std::ostream& out, std::string_view s) {
    out << '"';
    for (unsigned char c : s) {
        if (c == '\\' || c == '"') out << '\\' << c;
        else if (c < 32) out << "\\u00" << "0123456789abcdef"[c >> 4] << "0123456789abcdef"[c & 15];
        else out << c;
    }
    out << '"';
}
template<class T> struct IsVector : std::false_type {};
template<class T, class Allocator> struct IsVector<std::vector<T, Allocator>> : std::true_type {};
template<class Vector>
void list(std::ostream& out, const Vector& values) {
    out << '[';
    bool first = true;
    for (const auto& v : values) {
        if (!first) out << ',';
        first = false;
        if constexpr (IsVector<std::decay_t<decltype(v)>>::value) list(out, v);
        else if constexpr (std::is_floating_point_v<std::decay_t<decltype(v)>>) {
            if (v <= std::numeric_limits<double>::max() && v >= -std::numeric_limits<double>::max()) out << v;
            else out << "null";
        } else out << v;
    }
    out << ']';
}

struct Scope {
    static inline thread_local Scope* current = nullptr;
    Scope* parent = nullptr;
    Key key;
    std::uint64_t id = 0;
    std::int64_t start = 0, child = 0;
    int exceptions = 0;
    unsigned slot = 0;
    bool jacobian_recorded = false;
    std::array<Cell, 2> sweep{};
    const void* tape = nullptr;
    unsigned range = 0;
    explicit Scope(const char* name, unsigned first = 0, unsigned second = 0, unsigned third = 0)
        : key{name, stage, first, second, third, family, 0} {
        if (!active) return;
        parent = current;
        current = this;
        slot = trial;
        if (tracing) id = ++sequence;
        exceptions = std::uncaught_exceptions();
        start = now();
    }
    Scope(const Scope&) = delete;
    ~Scope() {
        if (!start) return;
        const auto elapsed = now() - start;
        if (std::uncaught_exceptions() > exceptions) key.outcome = -1;
        auto& cell = cells[key];
        ++cell.count;
        cell.inclusive += elapsed;
        cell.exclusive += elapsed - child;
        if (parent) parent->child += elapsed;
        current = parent;
        if (!tracing) return;
        // High-volume sweeps are retained as cells, and their numeric requests as events.
        const std::string_view name(key.name);
        if (name == "cppad_forward" || name == "cppad_reverse") {
            for (auto* outer = parent; outer; outer = outer->parent) {
                if (std::string_view(outer->key.name) != "ipopt_eval_h") continue;
                auto& sum = outer->sweep[name == "cppad_reverse"];
                ++sum.count;
                sum.exclusive += elapsed - child;
                // Retain sweep identity for the curvature call, where the read breakdown is priced.
                event("sweep", [&](auto& out) {
                    out << ',' << outer->id << ',' << id << ',' << (name == "cppad_reverse")
                        << ',' << key.first << ',' << key.second << ',' << key.third << ',' << range
                        << ',' << reinterpret_cast<std::uintptr_t>(tape) << ',' << elapsed-child;
                });
                break;
            }
            return;
        }
        if (name == "ipopt_eval_h" && !key.first) event("hessian_partition", [&](auto& out) {
            out << ',' << id << ',' << sweep[0].count << ',' << sweep[0].exclusive << ','
                << sweep[1].count << ',' << sweep[1].exclusive << ','
                << elapsed - sweep[0].exclusive - sweep[1].exclusive;
        });
        event("span", [&](auto& out) {
            out << ','; quoted(out, name);
            out << ',' << id << ',' << (parent ? parent->id : 0) << ',' << slot
                << ',' << start << ',' << elapsed << ',' << elapsed - child
                << ',' << key.stage << ",\"" << key.family << "\","
                << key.first << ',' << key.second << ',' << key.third << ',' << key.outcome;
        });
    }
    bool ok(bool result) {
        if (!result) key.outcome = -1;
        return result;
    }
    template<class Write>
    static void event(const char* type, Write write) {
        if (!tracing) return;
        std::ostringstream out;
        out << std::setprecision(17) << "[";
        quoted(out, type);
        write(out);
        out << ']';
        events.push_back(std::move(out).str());
    }
};
struct Family {
    char previous = family;
    Scope scope;
    static const char* label(char f) {
        switch (f) {
        case 'P': return "parent";
        case 'O': return "optimized";
        case 'S': return "supplied";
        default: return "split";
        }
    }
    explicit Family(char f) : scope((family = f, label(f))) {}
    ~Family() { family = previous; }
};
struct Trial {
    unsigned previous = trial;
    Scope scope;
    explicit Trial(bool supplied, bool retry)
        : scope((trial = active ? ++trial_sequence : 0, "start"), supplied, retry) {}
    ~Trial() { trial = previous; }
};

template<class Vector>
void state_key(const char* kind, std::initializer_list<double> scalar, const Vector& values,
               const void* model) {
    Scope::event(kind, [&](auto& out) {
        out << ',' << Scope::current->id << ',' << reinterpret_cast<std::uintptr_t>(model) << ',';
        list(out, scalar); out << ','; list(out, values);
    });
}
} // namespace ep235
