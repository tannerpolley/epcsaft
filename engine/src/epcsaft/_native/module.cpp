#include "eos.hpp"
#include "native_sdk.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <limits>
#include <utility>
#include <vector>

namespace py = pybind11;

PYBIND11_MODULE(_core, module) {
    py::class_<epcsaft::NativeDensityDiagnostics>(module, "NativeDensityDiagnostics")
        .def_readonly(
            "bracket_lower", &epcsaft::NativeDensityDiagnostics::bracket_lower
        )
        .def_readonly(
            "bracket_upper", &epcsaft::NativeDensityDiagnostics::bracket_upper
        )
        .def_readonly(
            "pressure_residual",
            &epcsaft::NativeDensityDiagnostics::pressure_residual
        )
        .def_readonly(
            "pressure_density_derivative",
            &epcsaft::NativeDensityDiagnostics::pressure_density_derivative
        )
        .def_readonly(
            "iterations", &epcsaft::NativeDensityDiagnostics::iterations
        )
        .def_readonly("stable", &epcsaft::NativeDensityDiagnostics::stable)
        .def_readonly("branch", &epcsaft::NativeDensityDiagnostics::branch);

    py::class_<epcsaft::NativeResult>(module, "NativeResult")
        .def_readonly("hard_chain", &epcsaft::NativeResult::hard_chain)
        .def_readonly("dispersion", &epcsaft::NativeResult::dispersion)
        .def_readonly("association", &epcsaft::NativeResult::association)
        .def_readonly("debye_huckel", &epcsaft::NativeResult::debye_huckel)
        .def_readonly("born_ssm_ds", &epcsaft::NativeResult::born_ssm_ds)
        .def_readonly(
            "residual_helmholtz", &epcsaft::NativeResult::residual_helmholtz
        )
        .def_readonly(
            "compressibility_factor",
            &epcsaft::NativeResult::compressibility_factor
        )
        .def_readonly("pressure", &epcsaft::NativeResult::pressure)
        .def_readonly(
            "residual_chemical_potential_over_rt",
            &epcsaft::NativeResult::residual_chemical_potential_over_rt
        )
        .def_readonly(
            "log_fugacity_coefficient",
            &epcsaft::NativeResult::log_fugacity_coefficient
        )
        .def_readonly(
            "residual_entropy", &epcsaft::NativeResult::residual_entropy
        )
        .def_readonly(
            "residual_internal_energy",
            &epcsaft::NativeResult::residual_internal_energy
        )
        .def_readonly(
            "residual_enthalpy", &epcsaft::NativeResult::residual_enthalpy
        )
        .def_readonly(
            "residual_isochoric_heat_capacity",
            &epcsaft::NativeResult::residual_isochoric_heat_capacity
        )
        .def_readonly(
            "residual_isobaric_heat_capacity",
            &epcsaft::NativeResult::residual_isobaric_heat_capacity
        )
        .def_readonly("molar_density", &epcsaft::NativeResult::molar_density)
        .def_readonly(
            "density_diagnostics", &epcsaft::NativeResult::density_diagnostics
        );

    py::class_<epcsaft::NativeReferenceSequence>(module, "NativeReferenceSequence")
        .def_readonly(
            "log_fugacity", &epcsaft::NativeReferenceSequence::log_fugacity
        )
        .def_readonly(
            "reference_molality",
            &epcsaft::NativeReferenceSequence::reference_molality
        )
        .def_readonly(
            "convergence_error",
            &epcsaft::NativeReferenceSequence::convergence_error
        );

    py::class_<epcsaft::NativeModel>(module, "NativeModel")
        .def(py::init([](
            std::vector<double> segment_count,
            std::vector<double> sigma_constant,
            std::vector<int> sigma_term_component,
            std::vector<double> sigma_term_amplitude,
            std::vector<double> sigma_term_exponent,
            std::vector<double> dispersion_energy_over_k,
            std::vector<double> k_ij,
            std::vector<double> l_ij,
            std::vector<int> site_component,
            std::vector<int> site_multiplicity,
            std::vector<int> association_active,
            std::vector<double> association_energy_over_k,
            std::vector<double> association_volume,
            std::vector<double> molar_mass_kg_per_mol,
            std::vector<int> charge_number,
            std::vector<double> relative_permittivity,
            double dielectric_ion_suppression_coefficient,
            std::vector<double> born_diameter_angstrom,
            std::vector<double> solvation_factor,
            double ionic_region_relative_permittivity,
            int electrolyte_model
        ) {
            return epcsaft::NativeModel(epcsaft::NativeParameters{
                std::move(segment_count),
                std::move(sigma_constant),
                std::move(sigma_term_component),
                std::move(sigma_term_amplitude),
                std::move(sigma_term_exponent),
                std::move(dispersion_energy_over_k),
                std::move(k_ij),
                std::move(l_ij),
                std::move(site_component),
                std::move(site_multiplicity),
                std::move(association_active),
                std::move(association_energy_over_k),
                std::move(association_volume),
                std::move(molar_mass_kg_per_mol),
                std::move(charge_number),
                std::move(relative_permittivity),
                dielectric_ion_suppression_coefficient,
                std::move(born_diameter_angstrom),
                std::move(solvation_factor),
                ionic_region_relative_permittivity,
                electrolyte_model,
            });
        }))
        .def("evaluate", &epcsaft::NativeModel::evaluate)
        .def("evaluate_pressure", &epcsaft::NativeModel::evaluate_pressure)
        .def("reference_composition", &epcsaft::NativeModel::reference_composition)
        .def(
            "electrolyte_reference_sequence",
            [](const epcsaft::NativeModel& model,
               double temperature,
               double pressure,
               const std::vector<double>& solvent_fractions,
               const std::array<std::size_t, 2>& ion_indices,
               const std::array<int, 2>& stoichiometry,
               const std::string& phase) {
                return model.electrolyte_reference_sequence(
                    temperature,
                    pressure,
                    solvent_fractions,
                    ion_indices,
                    stoichiometry,
                    phase
                );
            }
        );

    module.def(
        "native_sdk_v1",
        &epcsaft::make_native_sdk_v1,
        py::arg("model"),
        py::arg("parameter_fingerprint"),
        py::arg("temperature_min_k"),
        py::arg("temperature_max_k"),
        py::arg("component_ids"),
        py::arg("ion_mole_fraction_max"),
        py::arg("phase_callbacks_admitted") = true,
        py::arg("active_born_admitted") = false,
        py::arg("aqueous_kij_admitted") = false,
        py::arg("pressure_min_pa") = std::numeric_limits<double>::quiet_NaN(),
        py::arg("pressure_max_pa") = std::numeric_limits<double>::quiet_NaN(),
        py::arg("neutral_reference_admitted") = false,
        py::arg("unknown_applicability_admitted") = false,
        py::arg("aqueous_solvation_factor_admitted") = false,
        py::arg("staged_trial_bundle_admitted") = false
    );
}
