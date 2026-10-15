#pragma once
#include "core/numerical_policy.hpp"
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

namespace epcsaft::bindings {
inline void solver_option(policy::SolverOptions& options, const std::string& name,
    nanobind::handle declared, unsigned scope = 63, bool zero_unset = false) {
    if (nanobind::isinstance<nanobind::bool_>(declared)
        || !(nanobind::isinstance<nanobind::float_>(declared) || nanobind::isinstance<nanobind::int_>(declared)))
        throw std::invalid_argument(name + " requires a finite positive number");
    const double value = nanobind::cast<double>(declared);
    if (!zero_unset || value != 0.0) options.set(name, value, scope);
}
inline policy::SolverOptions solver_options(const nanobind::object& mapping, unsigned scope = 63) {
    policy::SolverOptions options;
    if (!mapping.is_none()) {
        const auto entries = nanobind::cast<nanobind::dict>(nanobind::module_::import_("builtins").attr("dict")(mapping));
        for (const auto& item : entries)
            solver_option(options, nanobind::cast<std::string>(item.first), item.second, scope);
    }
    return options;
}
inline nanobind::list solver_options_view(const policy::SolverOptions& options, unsigned scope = 63) {
    nanobind::list records;
    for (std::size_t i = 0; i < options.values.size(); ++i) if (scope & (1u << i)) {
        nanobind::dict record;
        record["name"] = policy::SOLVER_POLICIES[i]->name; record["value"] = options.value(i);
        record["unit"] = policy::SOLVER_UNITS[i]; record["default"] = policy::SOLVER_POLICIES[i]->value;
        records.append(record);
    }
    return records;
}
}
#define EPCSAFT_SOLVER_OPTIONS(T, scope) .def_prop_ro("solver_options", [](const T& v) { return epcsaft::bindings::solver_options_view(v.solver_options, scope); })
