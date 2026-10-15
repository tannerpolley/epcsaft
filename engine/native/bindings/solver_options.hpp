#pragma once
#include "core/numerical_policy.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace epcsaft::bindings {
inline void solver_option(policy::SolverOptions& options, const std::string& name,
    pybind11::handle declared, unsigned scope = 63, bool zero_unset = false) {
    if (pybind11::isinstance<pybind11::bool_>(declared)
        || !(pybind11::isinstance<pybind11::float_>(declared) || pybind11::isinstance<pybind11::int_>(declared)))
        throw std::invalid_argument(name + " requires a finite positive number");
    const double value = pybind11::cast<double>(declared);
    if (!zero_unset || value != 0.0) options.set(name, value, scope);
}
inline policy::SolverOptions solver_options(const pybind11::object& mapping, unsigned scope = 63) {
    policy::SolverOptions options;
    if (!mapping.is_none()) for (const auto& item : pybind11::dict(mapping))
        solver_option(options, pybind11::cast<std::string>(item.first), item.second, scope);
    return options;
}
inline pybind11::list solver_options_view(const policy::SolverOptions& options, unsigned scope = 63) {
    pybind11::list records;
    for (std::size_t i = 0; i < options.values.size(); ++i) if (scope & (1u << i)) {
        pybind11::dict record;
        record["name"] = policy::SOLVER_POLICIES[i]->name; record["value"] = options.value(i);
        record["unit"] = policy::SOLVER_UNITS[i]; record["default"] = policy::SOLVER_POLICIES[i]->value;
        records.append(record);
    }
    return records;
}
}
#define EPCSAFT_SOLVER_OPTIONS(T, scope) .def_property_readonly("solver_options", [](const T& v) { return epcsaft::bindings::solver_options_view(v.solver_options, scope); })
