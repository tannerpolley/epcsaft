#pragma once

#include "surface/parameter_builder.hpp"
#include "properties/helmholtz_model.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

#include <utility>

namespace epcsaft::bindings {

template<class Result, class Check>
nanobind::object domain_result(Result result, const HelmholtzModel& model, Check check) {
    if (!model.domains) return nanobind::cast(std::move(result));
    std::vector<OutsideDomain> violations;
    check(violations, result);
    auto value = nanobind::cast(std::move(result));
    if (!violations.empty()) {
        auto rows = nanobind::cast(std::move(violations));
        value.attr("_domain_violations") = rows;
        nanobind::module_::import_("epcsaft.parameters").attr("_warn_domains")(rows);
    }
    return value;
}

template<class Result>
auto domain_point(Result result, const HelmholtzModel& model, const std::vector<double>& x,
                  double t, std::optional<double> p = std::nullopt) {
    return domain_result(std::move(result), model, [&](auto& rows, const Result&) {
        outside_domains(rows, model.parameters(), *model.domains, x, t, p);
    });
}

}  // namespace epcsaft::bindings

#define EPCSAFT_DOMAIN_VIOLATIONS .def_prop_ro("domain_violations", [](nanobind::handle result) { \
    return nanobind::getattr(result, "_domain_violations", nanobind::none()); })
