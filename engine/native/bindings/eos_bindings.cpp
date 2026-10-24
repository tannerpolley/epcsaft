#include "properties/helmholtz_model.hpp"
#include "surface/parameter_builder.hpp"

#include "fields.hpp"
#include "solver_options.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>

#include <algorithm>
#include <array>
#include <optional>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace nb = nanobind;

namespace {

// Match scalar __bool__ conversion; containers' __len__ alone is insufficient.
bool interval_inclusive(nb::handle value) {
    const auto* number = Py_TYPE(value.ptr())->tp_as_number;
    const int truth = value.is_none() ? 0 : number && number->nb_bool ? PyObject_IsTrue(value.ptr()) : -1;
    if (truth == 0 || truth == 1) return truth == 1;
    PyErr_Clear();
    throw nb::type_error("expected a Boolean-convertible scalar");
}

template <class T, class... Values, class... Args>
auto bind_coefficients(nb::module_& module, const char* name, Args... args) {
    const auto convert = [](nb::handle value) {
        const std::string_view type = Py_TYPE(value.ptr())->tp_name;
        if (nb::isinstance<nb::str>(value) || nb::isinstance<nb::bytes>(value)
            || !(PySequence_Check(value.ptr()) || PyGen_Check(value.ptr()) || PyAnySet_Check(value.ptr())
                || type == "dict_keys" || type == "dict_values" || type == "dict_items" || type == "map" || type == "zip"))
            throw nb::type_error("expected a coefficient sequence or iterable");
        std::vector<double> values;
        for (nb::handle item : nb::iter(value)) {
            if (!nb::try_cast(item, values.emplace_back())) throw nb::type_error("invalid coefficient value");
        }
        decltype(T::coefficients) result{};
        if constexpr (std::is_same_v<decltype(result), std::vector<double>>) result = std::move(values);
        else if (values.size() == result.size()) std::copy(values.begin(), values.end(), result.begin());
        else throw nb::type_error("invalid coefficient count");
        return result;
    };
    return nb::class_<T>(module, name)
        .def("__init__", [convert](T* self, nb::handle values, Values... rest) { new (self) T{convert(values), rest...}; },
            nb::arg("coefficients"), args...)
        .def_prop_rw("coefficients", [](const T& v) { return v.coefficients; },
            [convert](T& v, nb::handle values) { v.coefficients = convert(values); });
}

// Physical input records: every field is required and transported unchanged; the native
// model copies them at construction and native admission owns their consistency.
void register_thermal_inputs(nb::module_& module) {
    using namespace epcsaft;
    module.attr("GAS_CONSTANT_J_PER_MOL_K") = GAS_CONSTANT_J_PER_MOL_K;
    bind_coefficients<IdealPolynomial, double, double, double>(module, "IdealPolynomial", nb::arg("reference_temperature"),
        nb::arg("reference_enthalpy"), nb::arg("reference_entropy"))
        EPCSAFT_READWRITE(IdealPolynomial, reference_temperature, reference_enthalpy, reference_entropy);
    nb::class_<IdealShomate>(module, "IdealShomate")
        .def(nb::init<std::array<double, 8>, double>(), nb::arg("coefficients"), nb::arg("formation_enthalpy"))
        EPCSAFT_READWRITE(IdealShomate, coefficients, formation_enthalpy);
    bind_coefficients<IdealNASA7>(module, "IdealNASA7");
    nb::class_<IdealNASA9>(module, "IdealNASA9")
        .def(nb::init<std::array<double, 9>>(),
             nb::arg("coefficients"))
        .def_rw("coefficients", &IdealNASA9::coefficients);
    nb::class_<IdealDIPPR107>(module, "IdealDIPPR107")
        .def(nb::init<std::array<double, 5>, double, double, double>(), nb::arg("coefficients"), nb::arg("reference_temperature"),
            nb::arg("reference_enthalpy"), nb::arg("reference_entropy"))
        EPCSAFT_READWRITE(IdealDIPPR107, coefficients, reference_temperature, reference_enthalpy, reference_entropy);
    nb::class_<IdealCorrelation>(module, "IdealCorrelation")
        .def(nb::init<decltype(IdealCorrelation::form), double>(), nb::arg("form"), nb::arg("source_pressure"))
        EPCSAFT_READWRITE(IdealCorrelation, form, source_pressure);
    auto interval = nb::class_<IdealInterval>(module, "IdealInterval");
    interval.def("__init__", [](IdealInterval* self, double lower, double upper, nb::handle lo, nb::handle hi, IdealCorrelation corr) {
            new (self) IdealInterval{lower, upper, interval_inclusive(lo), interval_inclusive(hi), std::move(corr)};
        }, nb::arg("lower"), nb::arg("upper"), nb::arg("lower_inclusive").none(), nb::arg("upper_inclusive").none(), nb::arg("correlation"))
        EPCSAFT_READWRITE(IdealInterval, lower, upper, correlation);
    for (const auto& [name, field] : {std::pair{"lower_inclusive", &IdealInterval::lower_inclusive}, std::pair{"upper_inclusive", &IdealInterval::upper_inclusive}})
        interval.def_prop_rw(name, [field](const IdealInterval& v) { return v.*field; },
            [field](IdealInterval& v, nb::handle value) { v.*field = interval_inclusive(value); }, nb::for_setter(nb::arg("value").none()));
    nb::class_<ThermochemistryRecord>(module, "ThermochemistryRecord")
        .def(nb::init<double, std::vector<std::vector<IdealInterval>>>(), nb::arg("standard_pressure"), nb::arg("component_intervals"))
        EPCSAFT_READWRITE(ThermochemistryRecord, standard_pressure, component_intervals);
    nb::class_<ReactionLogPolynomial>(module, "ReactionLogPolynomial")
        .def(nb::init<double, double, double, double, double, double, double, std::string, double, double>(), nb::arg("a"), nb::arg("b"), nb::arg("c"), nb::arg("d"),
            nb::arg("reference_temperature"), nb::arg("temperature_min"),
            nb::arg("temperature_max"), nb::arg("standard_state_id"), nb::arg("e") = 0.0, nb::arg("f") = 0.0)
        EPCSAFT_READWRITE(ReactionLogPolynomial, a, b, c, d, reference_temperature, temperature_min, temperature_max,
            standard_state_id, e, f);
    nb::enum_<ParameterFamily>(module, "ParameterFamily")
        EPCSAFT_VALUES(ParameterFamily, BornCDielectric, IonicRegionRelativePermittivity, BornDiameter,
            DebyeHuckelDiameter, PackingDiameter, IonSuppressionCoefficient, SegmentCount, SigmaConstant,
            DispersionEnergyOverK, Kij, Lij, AssociationEnergyOverK, AssociationVolume, KijReciprocalTemperatureSlope,
            DielectricVolume, DielectricTemperature, TemperatureLawCoefficient, CationPairMixing, CationTripleMixing);
    // Native slot indices of one coordinate; a reciprocal pair lists both of its slots.
    nb::class_<ActiveParameter>(module, "ActiveParameter")
        .def(nb::init<ParameterFamily, std::vector<std::size_t>>(), nb::arg("family"), nb::arg("indices"))
        EPCSAFT_READWRITE(ActiveParameter, family, indices);
}

}  // namespace

void register_eos(nb::module_& module) {
    register_thermal_inputs(module);
    nb::class_<epcsaft::ParameterEntry>(module, "ParameterEntry")
        .def(nb::init<>())
        EPCSAFT_READWRITE(epcsaft::ParameterEntry, identity, field, index, magnitude, unit);
    nb::class_<epcsaft::ParameterDocument>(module, "ParameterDocument")
        .def(nb::init<>())
        EPCSAFT_READWRITE(epcsaft::ParameterDocument, component_ids, site_component, site_multiplicity, base_choice,
            association_choice, electrolyte_choice, ionic_dispersion_choice, permittivity_choice, polar_choice,
            entries);
    module.def("build_parameters", &epcsaft::build_parameters);
    module.def("prepare_thermochemistry", &epcsaft::prepare_thermochemistry, nb::arg("record"), nb::arg("source_gas_constant"));
    module.def("fitted_document", &epcsaft::fitted_document);

    nb::class_<epcsaft::Parameters>(module, "Parameters")
        .def_prop_ro("component_count", &epcsaft::Parameters::component_count)
        .def_prop_ro("site_count", &epcsaft::Parameters::site_count)
        .def_ro("ideal_gas", &epcsaft::Parameters::ideal_gas)
        .def_ro("molar_masses", &epcsaft::Parameters::molar_mass_kg_per_mol, "kg/mol, in component_ids order")
        .def_ro("temperature_law_coefficient_identities", &epcsaft::Parameters::temperature_law_coefficient_identity,
            "source identity of each TemperatureLawCoefficient slot");
    nb::enum_<epcsaft::DerivativeOrder>(module, "DerivativeOrder")
        EPCSAFT_VALUES(epcsaft::DerivativeOrder, Value, Gradient, PressureSlope, Hessian);
    nb::enum_<epcsaft::NeutralProperty>(module, "NeutralProperty")
        EPCSAFT_VALUES(epcsaft::NeutralProperty, ResidualHelmholtz, Pressure);
    nb::class_<epcsaft::NeutralPropertyDerivatives>(module, "NeutralPropertyDerivatives")
        EPCSAFT_READONLY(epcsaft::NeutralPropertyDerivatives, value, gradient, hessian);
    nb::enum_<epcsaft::PropertyObservable>(module, "PropertyObservable")
        EPCSAFT_VALUES(epcsaft::PropertyObservable, IdealEnthalpy, IdealInternalEnergy, IdealEntropy,
            IdealIsobaricHeatCapacity, IdealIsochoricHeatCapacity, ResidualEnthalpy, ResidualInternalEnergy,
            ResidualEntropy, ResidualIsobaricHeatCapacity, ResidualIsochoricHeatCapacity, TotalEnthalpy,
            TotalInternalEnergy, TotalEntropy, TotalIsobaricHeatCapacity, TotalIsochoricHeatCapacity,
            LogFugacityCoefficient);
    nb::enum_<epcsaft::PropertyActionStatus>(module, "PropertyActionStatus")
        EPCSAFT_VALUES(epcsaft::PropertyActionStatus, Available, MissingThermochemistry, OutsideInterval,
            ThermalBoundary, SupportChange, IonicInsertion, AssociationSingular, AssociationIllConditioned,
            AssociationResidual, DensityUnavailable, PressureSingular, CentralPropertyUnavailable, Nonfinite);
    // A native family/slot coordinate increment in its ALG-PARAMETERS unit.
    nb::class_<epcsaft::ParameterActionDirection>(module, "ParameterActionDirection")
        .def(nb::init<epcsaft::ParameterFamily, std::vector<std::size_t>, double>(), nb::arg("family"), nb::arg("indices"), nb::arg("value"))
        EPCSAFT_READONLY(epcsaft::ParameterActionDirection, family, indices, value);
    nb::class_<epcsaft::DensityActionDirection>(module, "DensityActionDirection")
        .def(nb::init<double, double, std::vector<double>, std::vector<epcsaft::ParameterActionDirection>>(), nb::arg("temperature") = 0.0, nb::arg("molar_density") = 0.0,
            nb::arg("composition") = std::vector<double>{},
            nb::arg("parameter") = std::vector<epcsaft::ParameterActionDirection>{})
        EPCSAFT_READONLY(epcsaft::DensityActionDirection, temperature, molar_density, composition, parameter);
    nb::class_<epcsaft::PressureActionDirection>(module, "PressureActionDirection")
        .def(nb::init<double, double, std::vector<double>, std::vector<epcsaft::ParameterActionDirection>>(), nb::arg("temperature") = 0.0, nb::arg("pressure") = 0.0,
            nb::arg("composition") = std::vector<double>{},
            nb::arg("parameter") = std::vector<epcsaft::ParameterActionDirection>{})
        EPCSAFT_READONLY(epcsaft::PressureActionDirection, temperature, pressure, composition, parameter);
    nb::class_<epcsaft::PropertyActionResult>(module, "PropertyActionResult")
        EPCSAFT_READONLY(epcsaft::PropertyActionResult, observable, component, value, action, status, energy_order);
    nb::enum_<epcsaft::DensityBranch>(module, "DensityBranch")
        EPCSAFT_VALUES(epcsaft::DensityBranch, Unspecified, Vapor, Liquid)
        .def("__int__", [](epcsaft::DensityBranch value) { return static_cast<int>(value); });
    using Diagnostics = epcsaft::PropertyActionDiagnostics;
    nb::class_<Diagnostics>(module, "PropertyActionDiagnostics")
        EPCSAFT_READONLY(Diagnostics, site_rank, site_factorizations, site_solves, site_pivot_ratio,
            association_residual, pressure_density_derivative, pressure_residual, stability_determinant,
            pressure_factorizations, pressure_solves, thermal_evaluations, recordings, directional_paths,
            forward_sweeps, branch, stable_root_count);
    nb::class_<epcsaft::ResidualEnergyAction>(module, "ResidualEnergyAction")
        EPCSAFT_READONLY(epcsaft::ResidualEnergyAction, action, status, diagnostics);
    nb::class_<epcsaft::PropertyActionBatch>(module, "PropertyActionBatch")
        EPCSAFT_SOLVER_OPTIONS(epcsaft::PropertyActionBatch, 7)
        EPCSAFT_READONLY(epcsaft::PropertyActionBatch, results, maximum_energy_order, diagnostics, molar_density,
            density_action, density_mixed_action, density_message, density_convergence, density_status);
    nb::enum_<epcsaft::ThermochemistryStatus>(module, "ThermochemistryStatus")
        EPCSAFT_VALUES(epcsaft::ThermochemistryStatus, MissingRecord, Available, OutsideInterval, Nonfinite);
    nb::class_<epcsaft::ReactionThermochemistry>(module, "ReactionThermochemistry")
        EPCSAFT_READONLY(epcsaft::ReactionThermochemistry, gibbs, enthalpy, heat_capacity, pressure_gibbs, ln_k_pressure, ln_k_density,
            k_pressure, k_density, gibbs_derivatives, status, cause);
    module.def("reaction_thermochemistry_density", &epcsaft::reaction_thermochemistry_density,
        nb::arg("record"), nb::arg("stoichiometry"), nb::arg("temperature"));
    nb::class_<epcsaft::State>(module, "State", "One homogeneous state in SI molar units (K, Pa, mol/m^3, J/mol); component lists in component_ids order.")
        .def_ro("requested_density_branch", &epcsaft::State::requested_density_branch)
        .def_prop_ro("returned_density_branch", &epcsaft::State::returned_density_branch)
        .def_prop_ro("density_branch_mismatch", &epcsaft::State::density_branch_mismatch)
        EPCSAFT_SOLVER_OPTIONS(epcsaft::State, 7)
        EPCSAFT_READONLY(epcsaft::State, hard_chain, dispersion, association, debye_huckel, born, dipolar, cation_mixing,
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
            density_branch, stable_root_count)
        EPCSAFT_READONLY(epcsaft::State, isothermal_compressibility, isobaric_expansivity, heat_capacity_ratio,
            isentropic_compressibility, speed_of_sound, joule_thomson_coefficient, partial_molar_volume,
            partial_molar_enthalpy, partial_molar_entropy);
    nb::class_<epcsaft::DensityDiagnosticPoint>(module, "DensityDiagnosticPoint")
        EPCSAFT_READONLY(epcsaft::DensityDiagnosticPoint, molar_density, packing_fraction, pressure,
            pressure_density_derivative, lower, upper, residual, mechanical_sign, linked, declared, cause, stop);
    nb::class_<epcsaft::DensityDiagnosticBracket>(module, "DensityDiagnosticBracket")
        EPCSAFT_READONLY(epcsaft::DensityDiagnosticBracket, lower, upper, lower_residual, upper_residual);
    nb::class_<epcsaft::DensityDiagnostics>(module, "DensityDiagnostics")
        EPCSAFT_SOLVER_OPTIONS(epcsaft::DensityDiagnostics, 7)
        EPCSAFT_READONLY(epcsaft::DensityDiagnostics, samples, pressure_roots, spinodals, declared_brackets, unresolved,
            evaluation_order, finite_domain, geometric_cap, declared_crossings, evaluations, stop);
    nb::class_<epcsaft::SweepCounters>(module, "SweepCounters")
        EPCSAFT_READONLY(epcsaft::SweepCounters, forward, reverse, hessian_calls, recordings);
    nb::class_<epcsaft::HelmholtzModel>(module, "HelmholtzModel")
        .def("__init__", [](epcsaft::HelmholtzModel* self, epcsaft::Parameters p, std::vector<epcsaft::ActiveParameter> a,
                         std::optional<epcsaft::ThermochemistryRecord> t, nb::object options) {
            new (self) epcsaft::HelmholtzModel(std::move(p), std::move(a), std::move(t), epcsaft::bindings::solver_options(options, 7));
        },
             nb::arg("parameters"), nb::arg("active_parameters") = std::vector<epcsaft::ActiveParameter>{},
             nb::arg("thermochemistry") = std::nullopt, nb::arg("solver_options") = nb::none())
        .def_prop_ro("solver_options", [](const epcsaft::HelmholtzModel& v) { return epcsaft::bindings::solver_options_view(v.solver_options(), 7); })
        .def("state", [](const epcsaft::HelmholtzModel& model, double temperature, double molar_density,
                         const std::vector<double>& composition, bool derived_properties) {
            model.validate_electroneutrality(composition);
            return model.state(temperature, molar_density, composition, derived_properties);
        }, nb::arg("temperature"), nb::arg("molar_density"), nb::arg("mole_fractions"),
             nb::arg("derived_properties") = false)
        .def("neutral_property", &epcsaft::HelmholtzModel::neutral_property)
        .def("density_actions", [](const epcsaft::HelmholtzModel& model,
                                    double temperature, double molar_density,
                                    const std::vector<double>& composition,
                                    const std::vector<epcsaft::PropertyObservable>& observables,
                                    const std::vector<epcsaft::DensityActionDirection>& directions) {
            epcsaft::core::Scope scope("density_actions");
            model.validate_electroneutrality(composition);
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
            model.validate_electroneutrality(composition);
            return model.pressure_actions(temperature, pressure, composition, branch,
                observables, directions, anchor_density);
        },
             nb::arg("temperature"), nb::arg("pressure"), nb::arg("mole_fractions"),
             nb::arg("branch"), nb::arg("observables"), nb::arg("directions"),
             nb::arg("anchor_density") = std::nullopt)
        .def("packing_fraction_coefficient", [](const epcsaft::HelmholtzModel& model,
                                               double temperature, const std::vector<double>& composition) {
            epcsaft::core::Scope scope("packing_fraction_coefficient");
            return model.packing_fraction_coefficient(temperature, composition);
        })
        .def("state_at_pressure", [](const epcsaft::HelmholtzModel& model, double temperature, double pressure,
                                     const std::vector<double>& composition, epcsaft::DensityBranch branch,
                                     std::optional<double> anchor_density, bool derived_properties) {
            model.validate_electroneutrality(composition);
            return model.state_at_pressure(temperature, pressure, composition, branch, anchor_density,
                derived_properties);
        },
             nb::arg("temperature"), nb::arg("pressure"), nb::arg("mole_fractions"),
             nb::arg("branch"), nb::arg("anchor_density") = std::nullopt, nb::arg("derived_properties") = false)
        .def("molar_density_upper_bound", &epcsaft::HelmholtzModel::molar_density_upper_bound)
        .def("density_diagnostics", &epcsaft::HelmholtzModel::density_diagnostics,
             nb::arg("temperature"), nb::arg("pressure"), nb::arg("composition"), nb::arg("nodes"),
             nb::arg("density_tolerance"), nb::arg("log_tolerance"), nb::arg("max_evaluations"),
             nb::arg("max_seconds"), nb::arg("adaptive"), nb::arg("finite_anchor"),
             nb::arg("finite_upper"), nb::arg("evaluation_order"))
        .def("counters", &epcsaft::HelmholtzModel::counters)
        .def("reset_counters", &epcsaft::HelmholtzModel::reset_counters);
    nb::class_<epcsaft::InfiniteDilutionReference>(module, "InfiniteDilutionReference")
        EPCSAFT_SOLVER_OPTIONS(epcsaft::InfiniteDilutionReference, 7)
        EPCSAFT_READONLY(epcsaft::InfiniteDilutionReference, log_fugacity, molality, convergence_error,
            parameter_derivative, pure_solvent_log_fugacity, pure_solvent_state, solvent_state, diagnostics,
            solvent_status, solvent_cause, solvent_convergence);
    nb::class_<epcsaft::SolventActivity>(module, "SolventActivity")
        EPCSAFT_READONLY(epcsaft::SolventActivity, solvent_activity, osmotic_coefficient,
            solvent_activity_parameter_derivative, osmotic_parameter_derivative, parameter_status, parameter_cause,
            status, cause, osmotic_cause, state);
    nb::class_<epcsaft::HenryConstant>(module, "HenryConstant")
        EPCSAFT_READONLY(epcsaft::HenryConstant, component_index, raw_log_reference, henry_constant_px,
            observed_terminal_change, status, message);
    module.def("henry_constants", &epcsaft::henry_constants);
    module.def("solvent_reference", [](const epcsaft::HelmholtzModel& model, double T, double P, const std::vector<double>& x) {
        return epcsaft::solvent_reference(model, T, P, x);
    });
    module.def("solvent_activity", [](const epcsaft::HelmholtzModel& model, double T, double P,
        const std::vector<double>& x, const epcsaft::InfiniteDilutionReference& reference) {
        return epcsaft::solvent_activity(model, T, P, x, reference);
    });
    nb::class_<epcsaft::ElectrolyteActivity>(module, "ElectrolyteActivity", "One salt at molality m (mol/kg solvent), infinite-dilution convention; mean_ionic_molality_scale is the molality-scale gamma+-.")
        EPCSAFT_SOLVER_OPTIONS(epcsaft::ElectrolyteActivity, 7)
        EPCSAFT_READONLY(epcsaft::ElectrolyteActivity, single_ion, mean_ionic_mole_fraction_scale,
            mean_ionic_log_molality_scale, mean_ionic_molality_scale, osmotic_coefficient, mean_ionic_parameter_derivative);
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
