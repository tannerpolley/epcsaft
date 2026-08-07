#include <array>
#include <cstddef>
#include <iostream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "native_transport.hpp"

#include "held2_certificates.hpp"
#include "flash.hpp"
#include "result_json.hpp"

namespace py = pybind11;

namespace epcsaft_engine::equilibrium {
void bind_chemical_equilibrium(py::module_& module);
}

namespace {

const epcsaft_native_transport& checked_transport(const py::capsule& capsule) {
    const char* name = capsule.name();
    if (name == nullptr
        || std::string_view(name) != EPCSAFT_NATIVE_TRANSPORT_CAPSULE_NAME) {
        throw py::value_error("expected private Engine native transport");
    }
    const auto* transport = capsule.get_pointer<epcsaft_native_transport>();
    if (transport == nullptr || transport->model_context == nullptr
        || transport->evaluate_pure_phase == nullptr) {
        throw py::value_error("native transport is incomplete");
    }
    return *transport;
}

const epcsaft_native_transport& checked_mixture_transport(
    const py::capsule& capsule
) {
    const epcsaft_native_transport& transport = checked_transport(capsule);
    if (transport.component_count != 2) {
        throw py::value_error("provider capsule mixture component count must be two");
    }
    if (transport.evaluate_mixture_phase == nullptr) {
        throw py::value_error("provider capsule is missing its mixture phase evaluator");
    }
    return transport;
}

py::dict transport_info(const py::capsule& capsule) {
    const epcsaft_native_transport& transport = checked_transport(capsule);
    py::dict result;
    result["capsule_name"] = EPCSAFT_NATIVE_TRANSPORT_CAPSULE_NAME;
    result["has_model_context"] = transport.model_context != nullptr;
    result["has_evaluate_pure_phase"] = transport.evaluate_pure_phase != nullptr;
    result["component_count"] = transport.component_count;
    result["has_evaluate_mixture_phase"] =
        transport.evaluate_mixture_phase != nullptr;
    return result;
}

py::dict mixture_phase_to_dict(
    const epcsaft_engine::equilibrium::MixturePhaseEvaluation& phase
) {
    py::dict result;
    result["helmholtz_over_rt_reference_amount"] = phase.value;
    result["gradient"] = phase.gradient;
    result["hessian"] = phase.hessian;
    result["pressure_pa"] = phase.pressure_pa;
    result["parameter_fingerprint"] = phase.parameter_fingerprint;
    return result;
}

py::dict evaluate_mixture_phase(
    const py::capsule& capsule,
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3,
    const std::string& expected_fingerprint
) {
    const epcsaft_engine::equilibrium::ProviderContext provider(
        checked_mixture_transport(capsule), expected_fingerprint
    );
    return mixture_phase_to_dict(
        provider.evaluate_mixture(temperature_k, amounts_mol, volume_m3)
    );
}

py::dict evaluate_phase(
    const py::capsule& capsule,
    double temperature_k,
    double amount_mol,
    double volume_m3,
    const std::string& expected_fingerprint
) {
    const epcsaft_native_transport& transport = checked_transport(capsule);
    const epcsaft_engine::equilibrium::ProviderContext provider(
        transport, expected_fingerprint
    );
    const epcsaft_engine::equilibrium::PhaseEvaluation evaluation =
        provider.evaluate(temperature_k, amount_mol, volume_m3);
    const epcsaft_phase_block_result& phase = evaluation.provider;

    py::dict result;
    result["status"] = phase.status;
    result["amount_mol"] = amount_mol;
    result["volume_m3"] = volume_m3;
    result["helmholtz_over_rt_reference_amount"] = phase.helmholtz_over_rt_reference_amount;
    result["gradient"] = std::array<double, 2>{phase.gradient[0], phase.gradient[1]};
    result["hessian"] = std::array<double, 4>{
        phase.hessian[0], phase.hessian[1], phase.hessian[2], phase.hessian[3]
    };
    result["third"] = std::array<double, 8>{
        phase.third[0], phase.third[1], phase.third[2], phase.third[3],
        phase.third[4], phase.third[5], phase.third[6], phase.third[7]
    };
    result["pressure_pa"] = phase.pressure_pa;
    result["chemical_potential_over_rt"] = phase.chemical_potential_over_rt;
    result["parameter_fingerprint"] = evaluation.parameter_fingerprint;
    return result;
}

py::dict solve_tp_flash(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& overall_mole_fractions,
    const std::string& expected_fingerprint,
    bool trace
) {
    const epcsaft_engine::equilibrium::ProviderContext provider(
            checked_transport(capsule),
        expected_fingerprint
    );
    epcsaft_engine::equilibrium::Held2TerminalProgress terminal_progress(std::cout);
    epcsaft_engine::equilibrium::FlashResult flash;
    {
        py::gil_scoped_release release;
        flash = epcsaft_engine::equilibrium::solve_tp_flash(
            provider,
            {
                temperature_k,
                pressure_pa,
                overall_mole_fractions,
            },
            trace ? &terminal_progress : nullptr
        );
    }
    return py::module_::import("json")
        .attr("loads")(
            epcsaft_engine::equilibrium::flash_result_to_json(flash)
        ).cast<py::dict>();
}

py::dict held2_audit_farkas_certificate(
    const std::vector<std::vector<double>>& matrix,
    const std::vector<double>& row_lower,
    const std::vector<double>& row_upper,
    const std::vector<double>& column_lower,
    const std::vector<double>& column_upper,
    const std::vector<double>& row_ray
) {
    const auto audit =
        epcsaft_engine::equilibrium::audit_held2_farkas_certificate(
            matrix,
            row_lower,
            row_upper,
            column_lower,
            column_upper,
            row_ray
        );
    py::dict result;
    result["reason"] = audit.reason;
    result["row_multipliers"] = audit.row_multipliers;
    result["normalization"] = audit.normalization;
    result["row_sign_violation_inf"] = audit.row_sign_violation_inf;
    result["dual_feasibility_violation_inf"] =
        audit.dual_feasibility_violation_inf;
    result["contradiction_margin"] = audit.contradiction_margin;
    result["contradiction_scale"] = audit.contradiction_scale;
    result["contradiction_threshold"] = audit.contradiction_threshold;
    result["solver_ray_recovered_without_presolve"] =
        audit.solver_ray_recovered_without_presolve;
    result["accepted"] = audit.accepted;
    return result;
}

py::dict held2_audit_step5_kkt(
    const std::vector<double>& variables,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    const std::vector<double>& objective_gradient,
    const std::vector<std::vector<double>>& constraints,
    const std::vector<double>& constraint_upper,
    const std::vector<double>& lower_bound_multipliers,
    const std::vector<double>& upper_bound_multipliers,
    const std::vector<double>& constraint_multipliers,
    double pullback_residual,
    double pullback_scale,
    double pressure_residual,
    bool same_major_iteration,
    bool step4_binding_valid,
    bool pressure_branch_valid
) {
    const auto audit = epcsaft_engine::equilibrium::audit_held2_step5_kkt(
        variables,
        lower,
        upper,
        objective_gradient,
        constraints,
        constraint_upper,
        lower_bound_multipliers,
        upper_bound_multipliers,
        constraint_multipliers,
        pullback_residual,
        pullback_scale,
        pressure_residual,
        same_major_iteration,
        step4_binding_valid,
        pressure_branch_valid
    );
    py::dict result;
    result["reason"] = audit.reason;
    result["primal_residual_inf"] = audit.primal_residual_inf;
    result["dual_sign_violation_inf"] = audit.dual_sign_violation_inf;
    result["pullback_residual_inf"] = audit.pullback_residual_inf;
    result["pullback_scale"] = audit.pullback_scale;
    result["pressure_residual"] = audit.pressure_residual;
    result["stationarity_residual_inf"] =
        audit.stationarity_residual_inf;
    result["complementarity_inf"] = audit.complementarity_inf;
    result["active_constraint_count"] = audit.active_constraint_count;
    result["active_jacobian_rank"] = audit.active_jacobian_rank;
    result["accepted"] = audit.accepted;
    return result;
}

}  // namespace

PYBIND11_MODULE(_equilibrium, module) {
    module.doc() = "Native local equilibrium formulation";
    epcsaft_engine::equilibrium::bind_chemical_equilibrium(module);
    module.def("transport_info", &transport_info, py::arg("capsule"));
    module.def(
        "evaluate_mixture_phase",
        &evaluate_mixture_phase,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("amounts_mol"),
        py::arg("volume_m3"),
        py::arg("expected_fingerprint")
    );
    module.def(
        "_solve_tp_flash",
        &solve_tp_flash,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("overall_mole_fractions"),
        py::arg("expected_fingerprint"),
        py::kw_only(),
        py::arg("trace") = false
    );
    module.def(
        "_held2_audit_farkas_certificate",
        &held2_audit_farkas_certificate,
        py::arg("matrix"),
        py::arg("row_lower"),
        py::arg("row_upper"),
        py::arg("column_lower"),
        py::arg("column_upper"),
        py::arg("row_ray")
    );
    module.def(
        "_held2_audit_step5_kkt",
        &held2_audit_step5_kkt,
        py::arg("variables"),
        py::arg("lower"),
        py::arg("upper"),
        py::arg("objective_gradient"),
        py::arg("constraints"),
        py::arg("constraint_upper"),
        py::arg("lower_bound_multipliers"),
        py::arg("upper_bound_multipliers"),
        py::arg("constraint_multipliers"),
        py::arg("pullback_residual"),
        py::arg("pullback_scale"),
        py::arg("pressure_residual"),
        py::arg("same_major_iteration"),
        py::arg("step4_binding_valid"),
        py::arg("pressure_branch_valid")
    );
    module.def(
        "evaluate_phase",
        &evaluate_phase,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("amount_mol"),
        py::arg("volume_m3"),
        py::arg("expected_fingerprint")
    );
}
