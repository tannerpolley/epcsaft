#include "eos.hpp"
#include "native_transport_pybind.hpp"
#include "fit.hpp"
#include "phase_block.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <limits>
#include <utility>
#include <vector>

namespace py = pybind11;

void bind_pure_saturation(py::module_& module);

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
        .def_readonly("branch", &epcsaft::NativeDensityDiagnostics::branch)
        .def_readonly(
            "stable_root_count",
            &epcsaft::NativeDensityDiagnostics::stable_root_count
        )
        .def_readonly(
            "selected_stable_root_index",
            &epcsaft::NativeDensityDiagnostics::selected_stable_root_index
        );

    py::class_<epcsaft::NativeResult>(module, "NativeResult")
        .def_readonly("hard_chain", &epcsaft::NativeResult::hard_chain)
        .def_readonly("dispersion", &epcsaft::NativeResult::dispersion)
        .def_readonly("association", &epcsaft::NativeResult::association)
        .def_readonly("ares_polar", &epcsaft::NativeResult::ares_polar)
        .def_readonly("debye_huckel", &epcsaft::NativeResult::debye_huckel)
        .def_readonly("born_ssm_ds", &epcsaft::NativeResult::born_ssm_ds)
        .def_readonly(
            "bulk_relative_permittivity",
            &epcsaft::NativeResult::bulk_relative_permittivity
        )
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
            "pure_solvent_log_fugacity_coefficient",
            &epcsaft::NativeReferenceSequence::pure_solvent_log_fugacity_coefficient
        )
        .def_readonly(
            "reference_molality",
            &epcsaft::NativeReferenceSequence::reference_molality
        )
        .def_readonly(
            "convergence_error",
            &epcsaft::NativeReferenceSequence::convergence_error
        );

    py::class_<epcsaft::NativeIonSolvationDerivativeBlock>(
        module, "NativeIonSolvationDerivativeBlock"
    )
        .def_readonly(
            "solvation_gibbs_j_per_mol",
            &epcsaft::NativeIonSolvationDerivativeBlock::solvation_gibbs_j_per_mol
        )
        .def_readonly(
            "parameter_derivatives_j_per_mol",
            &epcsaft::NativeIonSolvationDerivativeBlock::
                parameter_derivatives_j_per_mol
        )
        .def_readonly(
            "reference_molality_mol_per_kg",
            &epcsaft::NativeIonSolvationDerivativeBlock::
                reference_molality_mol_per_kg
        )
        .def_readonly(
            "reference_convergence_error",
            &epcsaft::NativeIonSolvationDerivativeBlock::
                reference_convergence_error
        );

    py::class_<epcsaft::NativeAqueousMiacSolvationFactor>(
        module, "NativeAqueousMiacSolvationFactor"
    )
        .def_readonly(
            "log_mean_ionic_activity_coefficient_molality",
            &epcsaft::NativeAqueousMiacSolvationFactor::
                log_mean_ionic_activity_coefficient_molality
        )
        .def_readonly(
            "derivative", &epcsaft::NativeAqueousMiacSolvationFactor::derivative
        );

    py::class_<epcsaft::NativeIonFractionSuppression>(
        module, "NativeIonFractionSuppression"
    )
        .def_readonly(
            "bulk_relative_permittivity",
            &epcsaft::NativeIonFractionSuppression::bulk_relative_permittivity
        )
        .def_readonly(
            "derivative", &epcsaft::NativeIonFractionSuppression::derivative
        );

    py::class_<epcsaft::NativeParameterRecord>(module, "NativeParameterRecord")
        .def(
            py::init<std::string, double>(),
            py::arg("identity"),
            py::arg("value")
        )
        .def_readonly("identity", &epcsaft::NativeParameterRecord::identity)
        .def_readonly("value", &epcsaft::NativeParameterRecord::value);

    py::class_<epcsaft::NativePairCoordinate>(module, "NativePairCoordinate")
        .def(
            py::init<
                std::string,
                std::string,
                std::size_t,
                std::size_t,
                std::string,
                std::string,
                std::vector<std::array<std::size_t, 2>>,
                std::vector<double>,
                double
            >(),
            py::arg("identity"),
            py::arg("family"),
            py::arg("component_a"),
            py::arg("component_b"),
            py::arg("association_role_a"),
            py::arg("association_role_b"),
            py::arg("association_edges"),
            py::arg("association_energy_scales"),
            py::arg("value")
        )
        .def_readonly("identity", &epcsaft::NativePairCoordinate::identity)
        .def_readonly("family", &epcsaft::NativePairCoordinate::family)
        .def_readonly("component_a", &epcsaft::NativePairCoordinate::component_a)
        .def_readonly("component_b", &epcsaft::NativePairCoordinate::component_b)
        .def_readonly(
            "association_role_a", &epcsaft::NativePairCoordinate::association_role_a
        )
        .def_readonly(
            "association_role_b", &epcsaft::NativePairCoordinate::association_role_b
        )
        .def_readonly(
            "association_edges", &epcsaft::NativePairCoordinate::association_edges
        )
        .def_readonly(
            "association_energy_scales",
            &epcsaft::NativePairCoordinate::association_energy_scales
        )
        .def_readonly("value", &epcsaft::NativePairCoordinate::value);

    py::class_<epcsaft::NativePurePhaseBlock>(module, "NativePurePhaseBlock")
        .def_readonly(
            "helmholtz_over_rt_reference_amount",
            &epcsaft::NativePurePhaseBlock::helmholtz_over_rt_reference_amount
        )
        .def_readonly("gradient", &epcsaft::NativePurePhaseBlock::gradient)
        .def_readonly("hessian", &epcsaft::NativePurePhaseBlock::hessian)
        .def_readonly("third", &epcsaft::NativePurePhaseBlock::third)
        .def_readonly("pressure_pa", &epcsaft::NativePurePhaseBlock::pressure_pa)
        .def_readonly(
            "chemical_potential_over_rt",
            &epcsaft::NativePurePhaseBlock::chemical_potential_over_rt
        );

    py::class_<epcsaft::NativeFixedTopologyPhaseBlock>(
        module, "NativeFixedTopologyPhaseBlock"
    )
        .def_readonly(
            "helmholtz_over_rt_reference_amount",
            &epcsaft::NativeFixedTopologyPhaseBlock::
                helmholtz_over_rt_reference_amount
        )
        .def_readonly(
            "gradient", &epcsaft::NativeFixedTopologyPhaseBlock::gradient
        )
        .def_readonly(
            "hessian", &epcsaft::NativeFixedTopologyPhaseBlock::hessian
        )
        .def_readonly(
            "pressure_pa", &epcsaft::NativeFixedTopologyPhaseBlock::pressure_pa
        )
        .def_readonly(
            "chemical_potential_over_rt",
            &epcsaft::NativeFixedTopologyPhaseBlock::chemical_potential_over_rt
        );

    py::class_<epcsaft::NativeMixturePhaseBlock>(module, "NativeMixturePhaseBlock")
        .def_readonly(
            "helmholtz_over_rt_reference_amount",
            &epcsaft::NativeMixturePhaseBlock::helmholtz_over_rt_reference_amount
        )
        .def_readonly("gradient", &epcsaft::NativeMixturePhaseBlock::gradient)
        .def_readonly("hessian", &epcsaft::NativeMixturePhaseBlock::hessian)
        .def_readonly("pressure_pa", &epcsaft::NativeMixturePhaseBlock::pressure_pa);

    py::class_<epcsaft::NativeFugacityDerivativeBlock>(
        module, "NativeFugacityDerivativeBlock"
    )
        .def_readonly(
            "log_fugacity_coefficient",
            &epcsaft::NativeFugacityDerivativeBlock::log_fugacity_coefficient
        )
        .def_readonly(
            "pressure_derivative_per_pa",
            &epcsaft::NativeFugacityDerivativeBlock::pressure_derivative_per_pa
        )
        .def_readonly(
            "parameter_derivatives",
            &epcsaft::NativeFugacityDerivativeBlock::parameter_derivatives
        )
        .def_readonly(
            "molar_density_mol_per_m3",
            &epcsaft::NativeFugacityDerivativeBlock::molar_density_mol_per_m3
        )
        .def_readonly(
            "root_residual_pa",
            &epcsaft::NativeFugacityDerivativeBlock::root_residual_pa
        )
        .def_readonly(
            "pressure_density_derivative_pa_m3_per_mol",
            &epcsaft::NativeFugacityDerivativeBlock::
                pressure_density_derivative_pa_m3_per_mol
        )
        .def_readonly(
            "density_condition_number",
            &epcsaft::NativeFugacityDerivativeBlock::density_condition_number
        )
        .def_readonly("branch", &epcsaft::NativeFugacityDerivativeBlock::branch)
        .def_readonly(
            "bracket_lower_mol_per_m3",
            &epcsaft::NativeFugacityDerivativeBlock::bracket_lower_mol_per_m3
        )
        .def_readonly(
            "bracket_upper_mol_per_m3",
            &epcsaft::NativeFugacityDerivativeBlock::bracket_upper_mol_per_m3
        )
        .def_readonly(
            "stable_root_count",
            &epcsaft::NativeFugacityDerivativeBlock::stable_root_count
        )
        .def_readonly(
            "selected_stable_root_index",
            &epcsaft::NativeFugacityDerivativeBlock::selected_stable_root_index
        );

    py::class_<epcsaft::NativePackingFractionBlock>(
        module, "NativePackingFractionBlock"
    )
        .def_readonly(
            "packing_fraction", &epcsaft::NativePackingFractionBlock::packing_fraction
        )
        .def_readonly("gradient", &epcsaft::NativePackingFractionBlock::gradient)
        .def_readonly("hessian", &epcsaft::NativePackingFractionBlock::hessian);

    py::class_<epcsaft::NativeInversePackingGeometryBlock>(
        module, "NativeInversePackingGeometryBlock"
    )
        .def_readonly(
            "volume_m3", &epcsaft::NativeInversePackingGeometryBlock::volume_m3
        )
        .def_readonly(
            "gradient", &epcsaft::NativeInversePackingGeometryBlock::gradient
        )
        .def_readonly(
            "hessian", &epcsaft::NativeInversePackingGeometryBlock::hessian
        );

    py::class_<epcsaft::NativeMolarVolumeBoundsParameterBlock>(
        module, "NativeMolarVolumeBoundsParameterBlock"
    )
        .def_readonly(
            "bounds_m3_per_mol",
            &epcsaft::NativeMolarVolumeBoundsParameterBlock::bounds_m3_per_mol
        )
        .def_readonly(
            "parameter_derivatives_m3_per_mol",
            &epcsaft::NativeMolarVolumeBoundsParameterBlock::
                parameter_derivatives_m3_per_mol
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
            double ion_fraction_suppression_coefficient,
            std::vector<double> born_diameter_angstrom,
            std::vector<double> born_diameter_sigma_scale,
            std::vector<double> solvation_factor,
            double ionic_region_relative_permittivity,
            int electrolyte_model,
            std::vector<double> dielectric_volume_m3_per_mol,
            std::vector<double> dielectric_temperature_k,
            std::vector<double> ion_specific_suppression_coefficient,
            std::vector<double> dipole_scaling_debye2,
            std::vector<double> polarizability_scaling_angstrom3,
            std::vector<double> correlation_integral_parameter,
            int relative_permittivity_formulation,
            std::vector<double> dipole_moment_debye,
            std::vector<double> quadrupole_moment_debye_angstrom,
            int polar_formulation
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
                ion_fraction_suppression_coefficient,
                std::move(born_diameter_angstrom),
                std::move(born_diameter_sigma_scale),
                std::move(solvation_factor),
                ionic_region_relative_permittivity,
                electrolyte_model,
                std::move(dielectric_volume_m3_per_mol),
                std::move(dielectric_temperature_k),
                std::move(ion_specific_suppression_coefficient),
                std::move(dipole_scaling_debye2),
                std::move(polarizability_scaling_angstrom3),
                std::move(correlation_integral_parameter),
                relative_permittivity_formulation,
                std::move(dipole_moment_debye),
                std::move(quadrupole_moment_debye_angstrom),
                polar_formulation,
            });
        }))
        .def("evaluate", &epcsaft::NativeModel::evaluate)
        .def("evaluate_pressure", &epcsaft::NativeModel::evaluate_pressure)
        .def(
            "with_parameter_identities",
            &epcsaft::NativeModel::with_parameter_identities,
            py::arg("component_ids"),
            py::arg("site_ids")
        )
        .def(
            "evaluate_fixed_topology_association_active_parameters",
            &epcsaft::NativeModel::
                evaluate_fixed_topology_association_active_parameters,
            py::arg("temperature"),
            py::arg("amount"),
            py::arg("volume"),
            py::arg("active_parameters")
        )
        .def(
            "evaluate_fixed_topology_associating_mixture_active_parameters",
            &epcsaft::NativeModel::
                evaluate_fixed_topology_associating_mixture_active_parameters,
            py::arg("temperature"),
            py::arg("amounts"),
            py::arg("volume"),
            py::arg("active_parameters")
        )
        .def(
            "evaluate_mixture_phase_pair_coordinates",
            &epcsaft::NativeModel::evaluate_mixture_phase_pair_coordinates,
            py::arg("temperature"),
            py::arg("amounts"),
            py::arg("volume"),
            py::arg("coordinates")
        )
        .def(
            "evaluate_fugacity_derivatives_active",
            &epcsaft::NativeModel::evaluate_fugacity_derivatives_active,
            py::arg("temperature"),
            py::arg("pressure"),
            py::arg("composition"),
            py::arg("phase"),
            py::arg("active_parameters")
        )
        .def(
            "evaluate_ion_solvation_derivatives",
            &epcsaft::NativeModel::evaluate_ion_solvation_derivatives,
            py::arg("temperature"),
            py::arg("pressure"),
            py::arg("solvent_fractions"),
            py::arg("ion_indices"),
            py::arg("stoichiometry"),
            py::arg("phase"),
            py::arg("active_parameters")
        )
        .def(
            "evaluate_aqueous_miac_solvation_factor",
            &epcsaft::NativeModel::evaluate_aqueous_miac_solvation_factor
        )
        .def(
            "evaluate_ion_fraction_suppression",
            &epcsaft::NativeModel::evaluate_ion_fraction_suppression
        )
        .def("molar_volume_bounds", &epcsaft::NativeModel::molar_volume_bounds)
        .def(
            "molar_density_upper_bound",
            &epcsaft::NativeModel::molar_density_upper_bound
        )
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
        "native_transport",
        &epcsaft::make_native_transport,
        py::arg("model"),
        py::arg("parameter_fingerprint"),
        py::arg("temperature_min_k"),
        py::arg("temperature_max_k"),
        py::arg("component_ids"),
        py::arg("ion_mole_fraction_max"),
        py::arg("pressure_min_pa") = std::numeric_limits<double>::quiet_NaN(),
        py::arg("pressure_max_pa") = std::numeric_limits<double>::quiet_NaN(),
        py::arg("topology_fingerprint") = "",
        py::arg("component_order_fingerprint") = "",
        py::arg("domain_fingerprint") = "",
        py::arg("artifact_fingerprint") = ""
    );
    module.def(
        "source_reference_transfer",
        &epcsaft::evaluate_source_reference_transfer,
        py::arg("capsule"),
        py::arg("expected_parameter_fingerprint"),
        py::arg("expected_topology_fingerprint"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("component_ids"),
        py::arg("source_solvent_mole_fractions"),
        py::arg("reference_state_id"),
        py::arg("activity_convention_id"),
        py::arg("reference_state_fingerprint"),
        py::arg("standard_molality_mol_per_kg"),
        py::arg("source_reference_pressure_pa"),
        py::arg("required_derivatives")
    );
    bind_pure_saturation(module);
    bind_parameter_fit(module);
}
