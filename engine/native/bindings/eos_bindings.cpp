#include "properties/helmholtz_model.hpp"
#include "surface/parameter_builder.hpp"

#include "fields.hpp"
#include "solver_options.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>
#include <optional>
#include <utility>

namespace py = pybind11;

namespace {

// Physical input records: every field is required and transported unchanged; the native
// model copies them at construction and native admission owns their consistency.
void register_thermal_inputs(py::module_& module) {
    using namespace epcsaft;
    module.attr("GAS_CONSTANT_J_PER_MOL_K") = GAS_CONSTANT_J_PER_MOL_K;
    py::class_<IdealPolynomial>(module, "IdealPolynomial")
        .def(py::init<std::vector<double>, double, double, double>(), py::arg("coefficients"), py::arg("reference_temperature"),
            py::arg("reference_enthalpy"), py::arg("reference_entropy"))
        EPCSAFT_READWRITE(IdealPolynomial, coefficients, reference_temperature, reference_enthalpy, reference_entropy);
    py::class_<IdealShomate>(module, "IdealShomate")
        .def(py::init<std::array<double, 8>, double>(), py::arg("coefficients"), py::arg("formation_enthalpy"))
        EPCSAFT_READWRITE(IdealShomate, coefficients, formation_enthalpy);
    py::class_<IdealNASA7>(module, "IdealNASA7")
        .def(py::init<std::array<double, 7>>(),
             py::arg("coefficients"))
        .def_readwrite("coefficients", &IdealNASA7::coefficients);
    py::class_<IdealNASA9>(module, "IdealNASA9")
        .def(py::init<std::array<double, 9>>(),
             py::arg("coefficients"))
        .def_readwrite("coefficients", &IdealNASA9::coefficients);
    py::class_<IdealDIPPR107>(module, "IdealDIPPR107")
        .def(py::init<std::array<double, 5>, double, double, double>(), py::arg("coefficients"), py::arg("reference_temperature"),
            py::arg("reference_enthalpy"), py::arg("reference_entropy"))
        EPCSAFT_READWRITE(IdealDIPPR107, coefficients, reference_temperature, reference_enthalpy, reference_entropy);
    py::class_<IdealCorrelation>(module, "IdealCorrelation")
        .def(py::init<decltype(IdealCorrelation::form), double>(), py::arg("form"), py::arg("source_pressure"))
        EPCSAFT_READWRITE(IdealCorrelation, form, source_pressure);
    py::class_<IdealInterval>(module, "IdealInterval")
        .def(py::init<double, double, bool, bool, IdealCorrelation>(), py::arg("lower"), py::arg("upper"), py::arg("lower_inclusive"),
            py::arg("upper_inclusive"), py::arg("correlation"))
        EPCSAFT_READWRITE(IdealInterval, lower, upper, lower_inclusive, upper_inclusive, correlation);
    py::class_<ThermochemistryRecord>(module, "ThermochemistryRecord")
        .def(py::init<double, std::vector<std::vector<IdealInterval>>>(), py::arg("standard_pressure"), py::arg("component_intervals"))
        EPCSAFT_READWRITE(ThermochemistryRecord, standard_pressure, component_intervals);
    py::class_<ReactionLogPolynomial>(module, "ReactionLogPolynomial")
        .def(py::init<double, double, double, double, double, double, double, std::string>(), py::arg("a"), py::arg("b"), py::arg("c"), py::arg("d"),
            py::arg("reference_temperature"), py::arg("temperature_min"),
            py::arg("temperature_max"), py::arg("standard_state_id"))
        EPCSAFT_READWRITE(ReactionLogPolynomial, a, b, c, d, reference_temperature, temperature_min, temperature_max,
            standard_state_id);
    py::enum_<ParameterFamily>(module, "ParameterFamily")
        EPCSAFT_VALUES(ParameterFamily, BornCDielectric, IonicRegionRelativePermittivity, BornDiameter,
            DebyeHuckelDiameter, PackingDiameter, IonSuppressionCoefficient, SegmentCount, SigmaConstant,
            DispersionEnergyOverK, Kij, Lij, AssociationEnergyOverK, AssociationVolume, KijReciprocalTemperatureSlope,
            DielectricVolume, DielectricTemperature);
    // Native slot indices of one coordinate; a reciprocal pair lists both of its slots.
    py::class_<ActiveParameter>(module, "ActiveParameter")
        .def(py::init<ParameterFamily, std::vector<std::size_t>>(), py::arg("family"), py::arg("indices"))
        EPCSAFT_READWRITE(ActiveParameter, family, indices);
}

}  // namespace

void register_eos(py::module_& module) {
    register_thermal_inputs(module);
    py::class_<epcsaft::ParameterEntry>(module, "ParameterEntry")
        .def(py::init<>())
        EPCSAFT_READWRITE(epcsaft::ParameterEntry, identity, field, index, magnitude, unit);
    py::class_<epcsaft::ParameterDocument>(module, "ParameterDocument")
        .def(py::init<>())
        EPCSAFT_READWRITE(epcsaft::ParameterDocument, component_ids, site_component, site_multiplicity, base_choice,
            association_choice, electrolyte_choice, ionic_dispersion_choice, permittivity_choice, polar_choice,
            entries);
    module.def("build_parameters", &epcsaft::build_parameters);

    py::class_<epcsaft::Parameters>(module, "Parameters")
        .def_property_readonly("component_count", &epcsaft::Parameters::component_count)
        .def_property_readonly("site_count", &epcsaft::Parameters::site_count)
        .def_readonly("molar_masses", &epcsaft::Parameters::molar_mass_kg_per_mol, "kg/mol, in component_ids order");
    py::enum_<epcsaft::DerivativeOrder>(module, "DerivativeOrder")
        EPCSAFT_VALUES(epcsaft::DerivativeOrder, Value, Gradient, PressureSlope, Hessian);
    py::enum_<epcsaft::NeutralProperty>(module, "NeutralProperty")
        EPCSAFT_VALUES(epcsaft::NeutralProperty, ResidualHelmholtz, Pressure);
    py::class_<epcsaft::NeutralPropertyDerivatives>(module, "NeutralPropertyDerivatives")
        EPCSAFT_READONLY(epcsaft::NeutralPropertyDerivatives, value, gradient, hessian);
    py::enum_<epcsaft::PropertyObservable>(module, "PropertyObservable")
        EPCSAFT_VALUES(epcsaft::PropertyObservable, IdealEnthalpy, IdealInternalEnergy, IdealEntropy,
            IdealIsobaricHeatCapacity, IdealIsochoricHeatCapacity, ResidualEnthalpy, ResidualInternalEnergy,
            ResidualEntropy, ResidualIsobaricHeatCapacity, ResidualIsochoricHeatCapacity, TotalEnthalpy,
            TotalInternalEnergy, TotalEntropy, TotalIsobaricHeatCapacity, TotalIsochoricHeatCapacity,
            LogFugacityCoefficient);
    py::enum_<epcsaft::PropertyActionStatus>(module, "PropertyActionStatus")
        EPCSAFT_VALUES(epcsaft::PropertyActionStatus, Available, MissingThermochemistry, OutsideInterval,
            ThermalBoundary, SupportChange, IonicInsertion, AssociationSingular, AssociationIllConditioned,
            AssociationResidual, DensityUnavailable, PressureSingular, CentralPropertyUnavailable, Nonfinite);
    // A native family/slot coordinate increment in its ALG-PARAMETERS unit.
    py::class_<epcsaft::ParameterActionDirection>(module, "ParameterActionDirection")
        .def(py::init<epcsaft::ParameterFamily, std::vector<std::size_t>, double>(), py::arg("family"), py::arg("indices"), py::arg("value"))
        EPCSAFT_READONLY(epcsaft::ParameterActionDirection, family, indices, value);
    py::class_<epcsaft::DensityActionDirection>(module, "DensityActionDirection")
        .def(py::init<double, double, std::vector<double>, std::vector<epcsaft::ParameterActionDirection>>(), py::arg("temperature") = 0.0, py::arg("molar_density") = 0.0,
            py::arg("composition") = std::vector<double>{},
            py::arg("parameter") = std::vector<epcsaft::ParameterActionDirection>{})
        EPCSAFT_READONLY(epcsaft::DensityActionDirection, temperature, molar_density, composition, parameter);
    py::class_<epcsaft::PressureActionDirection>(module, "PressureActionDirection")
        .def(py::init<double, double, std::vector<double>, std::vector<epcsaft::ParameterActionDirection>>(), py::arg("temperature") = 0.0, py::arg("pressure") = 0.0,
            py::arg("composition") = std::vector<double>{},
            py::arg("parameter") = std::vector<epcsaft::ParameterActionDirection>{})
        EPCSAFT_READONLY(epcsaft::PressureActionDirection, temperature, pressure, composition, parameter);
    py::class_<epcsaft::PropertyActionResult>(module, "PropertyActionResult")
        EPCSAFT_READONLY(epcsaft::PropertyActionResult, observable, component, value, action, status, energy_order);
    py::enum_<epcsaft::DensityBranch>(module, "DensityBranch")
        EPCSAFT_VALUES(epcsaft::DensityBranch, Unspecified, Vapor, Liquid);
    using Diagnostics = epcsaft::PropertyActionDiagnostics;
    py::class_<Diagnostics>(module, "PropertyActionDiagnostics")
        EPCSAFT_READONLY(Diagnostics, site_rank, site_factorizations, site_solves, site_pivot_ratio,
            association_residual, pressure_density_derivative, pressure_residual, stability_determinant,
            pressure_factorizations, pressure_solves, thermal_evaluations, recordings, directional_paths,
            forward_sweeps, branch, stable_root_count);
    py::class_<epcsaft::PropertyActionBatch>(module, "PropertyActionBatch")
        EPCSAFT_SOLVER_OPTIONS(epcsaft::PropertyActionBatch, 7)
        EPCSAFT_READONLY(epcsaft::PropertyActionBatch, results, maximum_energy_order, diagnostics, molar_density,
            density_action, density_mixed_action, density_message, density_convergence, density_status);
    py::enum_<epcsaft::ThermochemistryStatus>(module, "ThermochemistryStatus")
        EPCSAFT_VALUES(epcsaft::ThermochemistryStatus, MissingRecord, Available, OutsideInterval, Nonfinite);
    py::class_<epcsaft::State>(module, "State", "One homogeneous state in SI molar units (K, Pa, mol/m^3, J/mol); component lists in component_ids order.")
        EPCSAFT_SOLVER_OPTIONS(epcsaft::State, 7)
        EPCSAFT_READONLY(epcsaft::State, hard_chain, dispersion, association, debye_huckel, born,
            bulk_relative_permittivity, residual_helmholtz, compressibility_factor, pressure,
            // None for a charged component on an ion-free support; see component_property_status (IonicInsertion).
            residual_chemical_potential_over_rt, component_property_status,
            log_fugacity_coefficient, log_fugacity_parameter_derivative, residual_entropy, residual_internal_energy,
            residual_enthalpy, residual_isochoric_heat_capacity, residual_isobaric_heat_capacity, ideal_enthalpy,
            ideal_internal_energy, ideal_entropy, ideal_isobaric_heat_capacity, ideal_isochoric_heat_capacity,
            total_enthalpy, total_internal_energy, total_entropy, total_isobaric_heat_capacity,
            total_isochoric_heat_capacity, molar_density, pressure_density_derivative, site_fractions,
            density_convergence, thermochemistry_status,
            // Unspecified and absent in the T-rho chart: no pressure-root search was performed.
            density_branch, stable_root_count);
    py::class_<epcsaft::DensityDiagnosticPoint>(module, "DensityDiagnosticPoint")
        EPCSAFT_READONLY(epcsaft::DensityDiagnosticPoint, molar_density, packing_fraction, pressure,
            pressure_density_derivative, lower, upper, residual, mechanical_sign, linked, declared, cause, stop);
    py::class_<epcsaft::DensityDiagnosticBracket>(module, "DensityDiagnosticBracket")
        EPCSAFT_READONLY(epcsaft::DensityDiagnosticBracket, lower, upper, lower_residual, upper_residual);
    py::class_<epcsaft::DensityDiagnostics>(module, "DensityDiagnostics")
        EPCSAFT_SOLVER_OPTIONS(epcsaft::DensityDiagnostics, 7)
        EPCSAFT_READONLY(epcsaft::DensityDiagnostics, samples, pressure_roots, spinodals, declared_brackets, unresolved,
            evaluation_order, finite_domain, geometric_cap, declared_crossings, evaluations, stop);
    py::class_<epcsaft::SweepCounters>(module, "SweepCounters")
        EPCSAFT_READONLY(epcsaft::SweepCounters, forward, reverse, hessian_calls, recordings);
    py::class_<epcsaft::HelmholtzModel>(module, "HelmholtzModel")
        .def(py::init([](epcsaft::Parameters p, std::vector<epcsaft::ActiveParameter> a,
                         std::optional<epcsaft::ThermochemistryRecord> t, py::object options) {
            return std::make_unique<epcsaft::HelmholtzModel>(std::move(p), std::move(a), std::move(t), epcsaft::bindings::solver_options(options, 7));
        }),
             py::arg("parameters"), py::arg("active_parameters") = std::vector<epcsaft::ActiveParameter>{},
             py::arg("thermochemistry") = std::nullopt, py::arg("solver_options") = py::none())
        .def_property_readonly("solver_options", [](const epcsaft::HelmholtzModel& v) { return epcsaft::bindings::solver_options_view(v.solver_options(), 7); })
        .def("state", &epcsaft::HelmholtzModel::state)
        .def("neutral_property", &epcsaft::HelmholtzModel::neutral_property)
        .def("density_actions", [](const epcsaft::HelmholtzModel& model,
                                    double temperature, double molar_density,
                                    const std::vector<double>& composition,
                                    const std::vector<epcsaft::PropertyObservable>& observables,
                                    const std::vector<epcsaft::DensityActionDirection>& directions) {
            epcsaft::core::Scope scope("density_actions");
            return model.density_actions(
                temperature, molar_density, composition, observables, directions);
        })
        .def("pressure_actions", [](const epcsaft::HelmholtzModel& model,
                                     double temperature, double pressure,
                                     const std::vector<double>& composition,
                                     epcsaft::DensityBranch branch,
                                     const std::vector<epcsaft::PropertyObservable>& observables,
                                     const std::vector<epcsaft::PressureActionDirection>& directions,
                                     std::optional<double> anchor_density) {
            epcsaft::core::Scope scope("pressure_actions");
            return model.pressure_actions(temperature, pressure, composition, branch,
                observables, directions, anchor_density);
        },
             py::arg("temperature"), py::arg("pressure"), py::arg("mole_fractions"),
             py::arg("branch"), py::arg("observables"), py::arg("directions"),
             py::arg("anchor_density") = std::nullopt)
        .def("packing_fraction_coefficient", [](const epcsaft::HelmholtzModel& model,
                                               double temperature, const std::vector<double>& composition) {
            epcsaft::core::Scope scope("packing_fraction_coefficient");
            return model.packing_fraction_coefficient(temperature, composition);
        })
        .def("state_at_pressure", &epcsaft::HelmholtzModel::state_at_pressure,
             py::arg("temperature"), py::arg("pressure"), py::arg("mole_fractions"),
             py::arg("branch"), py::arg("anchor_density") = std::nullopt)
        .def("molar_density_upper_bound", &epcsaft::HelmholtzModel::molar_density_upper_bound)
        .def("density_diagnostics", &epcsaft::HelmholtzModel::density_diagnostics,
             py::arg("temperature"), py::arg("pressure"), py::arg("composition"), py::arg("nodes"),
             py::arg("density_tolerance"), py::arg("log_tolerance"), py::arg("max_evaluations"),
             py::arg("max_seconds"), py::arg("adaptive"), py::arg("finite_anchor"),
             py::arg("finite_upper"), py::arg("evaluation_order"))
        .def("counters", &epcsaft::HelmholtzModel::counters)
        .def("reset_counters", &epcsaft::HelmholtzModel::reset_counters);
    py::class_<epcsaft::InfiniteDilutionReference>(module, "InfiniteDilutionReference")
        EPCSAFT_SOLVER_OPTIONS(epcsaft::InfiniteDilutionReference, 7)
        EPCSAFT_READONLY(epcsaft::InfiniteDilutionReference, log_fugacity, molality, convergence_error,
            parameter_derivative, parameter_convergence_error, pure_solvent_log_fugacity);
    py::class_<epcsaft::ElectrolyteActivity>(module, "ElectrolyteActivity", "One salt at molality m (mol/kg solvent), infinite-dilution convention; mean_ionic_molality_scale is the molality-scale gamma+-.")
        EPCSAFT_SOLVER_OPTIONS(epcsaft::ElectrolyteActivity, 7)
        EPCSAFT_READONLY(epcsaft::ElectrolyteActivity, single_ion, mean_ionic_mole_fraction_scale,
            mean_ionic_log_molality_scale, mean_ionic_molality_scale, osmotic_coefficient);
    module.def("infinite_dilution_reference", [](const epcsaft::HelmholtzModel& model, double temperature, double pressure,
        const std::vector<double>& solvent, const std::array<std::size_t, 2>& ions, const std::array<int, 2>& stoichiometry) {
        return epcsaft::infinite_dilution_reference(model, temperature, pressure, solvent, ions, stoichiometry);
    });
    module.def("electrolyte_activity", [](const epcsaft::HelmholtzModel& model, double temperature, double pressure,
        const std::vector<double>& solvent, const std::array<std::size_t, 2>& ions, const std::array<int, 2>& stoichiometry,
        double molality, const epcsaft::InfiniteDilutionReference& reference) {
        return epcsaft::electrolyte_activity(model, temperature, pressure, solvent, ions, stoichiometry, molality, reference);
    });
}
