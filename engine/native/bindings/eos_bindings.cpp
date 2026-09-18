#include "properties/helmholtz_model.hpp"
#include "surface/parameter_builder.hpp"

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
        .def(py::init([](std::vector<double> coefficients, double reference_temperature,
                         double reference_enthalpy, double reference_entropy) {
            return IdealPolynomial{std::move(coefficients), reference_temperature,
                                   reference_enthalpy, reference_entropy};
        }), py::arg("coefficients"), py::arg("reference_temperature"),
            py::arg("reference_enthalpy"), py::arg("reference_entropy"))
        .def_readwrite("coefficients", &IdealPolynomial::coefficients)
        .def_readwrite("reference_temperature", &IdealPolynomial::reference_temperature)
        .def_readwrite("reference_enthalpy", &IdealPolynomial::reference_enthalpy)
        .def_readwrite("reference_entropy", &IdealPolynomial::reference_entropy);
    py::class_<IdealShomate>(module, "IdealShomate")
        .def(py::init([](std::array<double, 8> coefficients, double formation_enthalpy) {
            return IdealShomate{coefficients, formation_enthalpy};
        }), py::arg("coefficients"), py::arg("formation_enthalpy"))
        .def_readwrite("coefficients", &IdealShomate::coefficients)
        .def_readwrite("formation_enthalpy", &IdealShomate::formation_enthalpy);
    py::class_<IdealNASA7>(module, "IdealNASA7")
        .def(py::init([](std::array<double, 7> coefficients) { return IdealNASA7{coefficients}; }),
             py::arg("coefficients"))
        .def_readwrite("coefficients", &IdealNASA7::coefficients);
    py::class_<IdealNASA9>(module, "IdealNASA9")
        .def(py::init([](std::array<double, 9> coefficients) { return IdealNASA9{coefficients}; }),
             py::arg("coefficients"))
        .def_readwrite("coefficients", &IdealNASA9::coefficients);
    py::class_<IdealDIPPR107>(module, "IdealDIPPR107")
        .def(py::init([](std::array<double, 5> coefficients, double reference_temperature,
                         double reference_enthalpy, double reference_entropy) {
            return IdealDIPPR107{coefficients, reference_temperature,
                                 reference_enthalpy, reference_entropy};
        }), py::arg("coefficients"), py::arg("reference_temperature"),
            py::arg("reference_enthalpy"), py::arg("reference_entropy"))
        .def_readwrite("coefficients", &IdealDIPPR107::coefficients)
        .def_readwrite("reference_temperature", &IdealDIPPR107::reference_temperature)
        .def_readwrite("reference_enthalpy", &IdealDIPPR107::reference_enthalpy)
        .def_readwrite("reference_entropy", &IdealDIPPR107::reference_entropy);
    py::class_<IdealCorrelation>(module, "IdealCorrelation")
        .def(py::init([](decltype(IdealCorrelation::form) form, double source_pressure) {
            return IdealCorrelation{std::move(form), source_pressure};
        }), py::arg("form"), py::arg("source_pressure"))
        .def_readwrite("form", &IdealCorrelation::form)
        .def_readwrite("source_pressure", &IdealCorrelation::source_pressure);
    py::class_<IdealInterval>(module, "IdealInterval")
        .def(py::init([](double lower, double upper, bool lower_inclusive, bool upper_inclusive,
                         IdealCorrelation correlation) {
            return IdealInterval{lower, upper, lower_inclusive, upper_inclusive,
                                 std::move(correlation)};
        }), py::arg("lower"), py::arg("upper"), py::arg("lower_inclusive"),
            py::arg("upper_inclusive"), py::arg("correlation"))
        .def_readwrite("lower", &IdealInterval::lower)
        .def_readwrite("upper", &IdealInterval::upper)
        .def_readwrite("lower_inclusive", &IdealInterval::lower_inclusive)
        .def_readwrite("upper_inclusive", &IdealInterval::upper_inclusive)
        .def_readwrite("correlation", &IdealInterval::correlation);
    py::class_<ThermochemistryRecord>(module, "ThermochemistryRecord")
        .def(py::init([](double standard_pressure,
                         std::vector<std::vector<IdealInterval>> component_intervals) {
            return ThermochemistryRecord{standard_pressure, std::move(component_intervals)};
        }), py::arg("standard_pressure"), py::arg("component_intervals"))
        .def_readwrite("standard_pressure", &ThermochemistryRecord::standard_pressure)
        .def_readwrite("component_intervals", &ThermochemistryRecord::component_intervals);
    py::class_<ReactionLogPolynomial>(module, "ReactionLogPolynomial")
        .def(py::init([](double a, double b, double c, double d, double reference_temperature,
                         double temperature_min, double temperature_max,
                         std::string standard_state_id) {
            return ReactionLogPolynomial{a, b, c, d, reference_temperature, temperature_min,
                                         temperature_max, std::move(standard_state_id)};
        }), py::arg("a"), py::arg("b"), py::arg("c"), py::arg("d"),
            py::arg("reference_temperature"), py::arg("temperature_min"),
            py::arg("temperature_max"), py::arg("standard_state_id"))
        .def_readwrite("a", &ReactionLogPolynomial::a)
        .def_readwrite("b", &ReactionLogPolynomial::b)
        .def_readwrite("c", &ReactionLogPolynomial::c)
        .def_readwrite("d", &ReactionLogPolynomial::d)
        .def_readwrite("reference_temperature", &ReactionLogPolynomial::reference_temperature)
        .def_readwrite("temperature_min", &ReactionLogPolynomial::temperature_min)
        .def_readwrite("temperature_max", &ReactionLogPolynomial::temperature_max)
        .def_readwrite("standard_state_id", &ReactionLogPolynomial::standard_state_id);
    py::enum_<ParameterFamily>(module, "ParameterFamily")
        .value("SegmentCount", ParameterFamily::SegmentCount)
        .value("SigmaConstant", ParameterFamily::SigmaConstant)
        .value("DispersionEnergyOverK", ParameterFamily::DispersionEnergyOverK)
        .value("Kij", ParameterFamily::Kij)
        .value("Lij", ParameterFamily::Lij)
        .value("AssociationEnergyOverK", ParameterFamily::AssociationEnergyOverK)
        .value("AssociationVolume", ParameterFamily::AssociationVolume);
    // Native slot indices of one coordinate; a reciprocal pair lists both of its slots.
    py::class_<ActiveParameter>(module, "ActiveParameter")
        .def(py::init([](ParameterFamily family, std::vector<std::size_t> indices) {
            return ActiveParameter{family, std::move(indices)};
        }), py::arg("family"), py::arg("indices"))
        .def_readwrite("family", &ActiveParameter::family)
        .def_readwrite("indices", &ActiveParameter::indices);
}

}  // namespace

void register_eos(py::module_& module) {
    register_thermal_inputs(module);
    py::class_<epcsaft::ParameterEntry>(module, "ParameterEntry")
        .def(py::init<>())
        .def_readwrite("identity", &epcsaft::ParameterEntry::identity)
        .def_readwrite("field", &epcsaft::ParameterEntry::field)
        .def_readwrite("index", &epcsaft::ParameterEntry::index)
        .def_readwrite("magnitude", &epcsaft::ParameterEntry::magnitude)
        .def_readwrite("unit", &epcsaft::ParameterEntry::unit);
    py::class_<epcsaft::ParameterDocument>(module, "ParameterDocument")
        .def(py::init<>())
        .def_readwrite("component_ids", &epcsaft::ParameterDocument::component_ids)
        .def_readwrite("site_component", &epcsaft::ParameterDocument::site_component)
        .def_readwrite("site_multiplicity", &epcsaft::ParameterDocument::site_multiplicity)
        .def_readwrite("electrolyte_choice", &epcsaft::ParameterDocument::electrolyte_choice)
        .def_readwrite("ionic_dispersion_choice", &epcsaft::ParameterDocument::ionic_dispersion_choice)
        .def_readwrite("permittivity_choice", &epcsaft::ParameterDocument::permittivity_choice)
        .def_readwrite("entries", &epcsaft::ParameterDocument::entries);
    module.def("build_parameters", &epcsaft::build_parameters);

    py::class_<epcsaft::Parameters>(module, "Parameters")
        .def_property_readonly("component_count", &epcsaft::Parameters::component_count)
        .def_property_readonly("site_count", &epcsaft::Parameters::site_count);
    py::enum_<epcsaft::DerivativeOrder>(module, "DerivativeOrder")
        .value("Value", epcsaft::DerivativeOrder::Value)
        .value("Gradient", epcsaft::DerivativeOrder::Gradient)
        .value("PressureSlope", epcsaft::DerivativeOrder::PressureSlope)
        .value("Hessian", epcsaft::DerivativeOrder::Hessian);
    py::enum_<epcsaft::NeutralProperty>(module, "NeutralProperty")
        .value("ResidualHelmholtz", epcsaft::NeutralProperty::ResidualHelmholtz)
        .value("Pressure", epcsaft::NeutralProperty::Pressure);
    py::class_<epcsaft::NeutralPropertyDerivatives>(module, "NeutralPropertyDerivatives")
        .def_readonly("value", &epcsaft::NeutralPropertyDerivatives::value)
        .def_readonly("gradient", &epcsaft::NeutralPropertyDerivatives::gradient)
        .def_readonly("hessian", &epcsaft::NeutralPropertyDerivatives::hessian);
    py::enum_<epcsaft::PropertyObservable>(module, "PropertyObservable")
        .value("IdealEnthalpy", epcsaft::PropertyObservable::IdealEnthalpy)
        .value("IdealInternalEnergy", epcsaft::PropertyObservable::IdealInternalEnergy)
        .value("IdealEntropy", epcsaft::PropertyObservable::IdealEntropy)
        .value("IdealIsobaricHeatCapacity", epcsaft::PropertyObservable::IdealIsobaricHeatCapacity)
        .value("IdealIsochoricHeatCapacity", epcsaft::PropertyObservable::IdealIsochoricHeatCapacity)
        .value("ResidualEnthalpy", epcsaft::PropertyObservable::ResidualEnthalpy)
        .value("ResidualInternalEnergy", epcsaft::PropertyObservable::ResidualInternalEnergy)
        .value("ResidualEntropy", epcsaft::PropertyObservable::ResidualEntropy)
        .value("ResidualIsobaricHeatCapacity", epcsaft::PropertyObservable::ResidualIsobaricHeatCapacity)
        .value("ResidualIsochoricHeatCapacity", epcsaft::PropertyObservable::ResidualIsochoricHeatCapacity)
        .value("TotalEnthalpy", epcsaft::PropertyObservable::TotalEnthalpy)
        .value("TotalInternalEnergy", epcsaft::PropertyObservable::TotalInternalEnergy)
        .value("TotalEntropy", epcsaft::PropertyObservable::TotalEntropy)
        .value("TotalIsobaricHeatCapacity", epcsaft::PropertyObservable::TotalIsobaricHeatCapacity)
        .value("TotalIsochoricHeatCapacity", epcsaft::PropertyObservable::TotalIsochoricHeatCapacity);
    py::enum_<epcsaft::PropertyActionStatus>(module, "PropertyActionStatus")
        .value("Available", epcsaft::PropertyActionStatus::Available)
        .value("MissingThermochemistry", epcsaft::PropertyActionStatus::MissingThermochemistry)
        .value("OutsideInterval", epcsaft::PropertyActionStatus::OutsideInterval)
        .value("ThermalBoundary", epcsaft::PropertyActionStatus::ThermalBoundary)
        .value("SupportChange", epcsaft::PropertyActionStatus::SupportChange)
        .value("IonicInsertion", epcsaft::PropertyActionStatus::IonicInsertion)
        .value("AssociationSingular", epcsaft::PropertyActionStatus::AssociationSingular)
        .value("AssociationIllConditioned", epcsaft::PropertyActionStatus::AssociationIllConditioned)
        .value("AssociationResidual", epcsaft::PropertyActionStatus::AssociationResidual)
        .value("DensityUnavailable", epcsaft::PropertyActionStatus::DensityUnavailable)
        .value("PressureSingular", epcsaft::PropertyActionStatus::PressureSingular)
        .value("CentralPropertyUnavailable", epcsaft::PropertyActionStatus::CentralPropertyUnavailable)
        .value("Nonfinite", epcsaft::PropertyActionStatus::Nonfinite);
    // A native family/slot coordinate increment in its ALG-PARAMETERS unit.
    py::class_<epcsaft::ParameterActionDirection>(module, "ParameterActionDirection")
        .def(py::init([](epcsaft::ParameterFamily family, std::vector<std::size_t> indices, double value) {
            return epcsaft::ParameterActionDirection{family, std::move(indices), value};
        }), py::arg("family"), py::arg("indices"), py::arg("value"))
        .def_readonly("family", &epcsaft::ParameterActionDirection::family)
        .def_readonly("indices", &epcsaft::ParameterActionDirection::indices)
        .def_readonly("value", &epcsaft::ParameterActionDirection::value);
    py::class_<epcsaft::DensityActionDirection>(module, "DensityActionDirection")
        .def(py::init([](double temperature, double molar_density, std::vector<double> composition,
                         std::vector<epcsaft::ParameterActionDirection> parameter) {
            return epcsaft::DensityActionDirection{
                temperature, molar_density, std::move(composition), std::move(parameter)};
        }), py::arg("temperature") = 0.0, py::arg("molar_density") = 0.0,
            py::arg("composition") = std::vector<double>{},
            py::arg("parameter") = std::vector<epcsaft::ParameterActionDirection>{})
        .def_readonly("temperature", &epcsaft::DensityActionDirection::temperature)
        .def_readonly("molar_density", &epcsaft::DensityActionDirection::molar_density)
        .def_readonly("composition", &epcsaft::DensityActionDirection::composition)
        .def_readonly("parameter", &epcsaft::DensityActionDirection::parameter);
    py::class_<epcsaft::PressureActionDirection>(module, "PressureActionDirection")
        .def(py::init([](double temperature, double pressure, std::vector<double> composition,
                         std::vector<epcsaft::ParameterActionDirection> parameter) {
            return epcsaft::PressureActionDirection{
                temperature, pressure, std::move(composition), std::move(parameter)};
        }), py::arg("temperature") = 0.0, py::arg("pressure") = 0.0,
            py::arg("composition") = std::vector<double>{},
            py::arg("parameter") = std::vector<epcsaft::ParameterActionDirection>{})
        .def_readonly("temperature", &epcsaft::PressureActionDirection::temperature)
        .def_readonly("pressure", &epcsaft::PressureActionDirection::pressure)
        .def_readonly("composition", &epcsaft::PressureActionDirection::composition)
        .def_readonly("parameter", &epcsaft::PressureActionDirection::parameter);
    py::class_<epcsaft::PropertyActionResult>(module, "PropertyActionResult")
        .def_readonly("observable", &epcsaft::PropertyActionResult::observable)
        .def_readonly("value", &epcsaft::PropertyActionResult::value)
        .def_readonly("action", &epcsaft::PropertyActionResult::action)
        .def_readonly("status", &epcsaft::PropertyActionResult::status)
        .def_readonly("energy_order", &epcsaft::PropertyActionResult::energy_order);
    py::enum_<epcsaft::DensityBranch>(module, "DensityBranch")
        .value("Unspecified", epcsaft::DensityBranch::Unspecified)
        .value("Vapor", epcsaft::DensityBranch::Vapor)
        .value("Liquid", epcsaft::DensityBranch::Liquid);
    using Diagnostics = epcsaft::PropertyActionDiagnostics;
    py::class_<Diagnostics>(module, "PropertyActionDiagnostics")
        .def_readonly("site_rank", &Diagnostics::site_rank)
        .def_readonly("site_factorizations", &Diagnostics::site_factorizations)
        .def_readonly("site_solves", &Diagnostics::site_solves)
        .def_readonly("site_pivot_ratio", &Diagnostics::site_pivot_ratio)
        .def_readonly("association_residual", &Diagnostics::association_residual)
        .def_readonly("pressure_density_derivative", &Diagnostics::pressure_density_derivative)
        .def_readonly("pressure_residual", &Diagnostics::pressure_residual)
        .def_readonly("stability_determinant", &Diagnostics::stability_determinant)
        .def_readonly("pressure_factorizations", &Diagnostics::pressure_factorizations)
        .def_readonly("pressure_solves", &Diagnostics::pressure_solves)
        .def_readonly("thermal_evaluations", &Diagnostics::thermal_evaluations)
        .def_readonly("recordings", &Diagnostics::recordings)
        .def_readonly("directional_paths", &Diagnostics::directional_paths)
        .def_readonly("forward_sweeps", &Diagnostics::forward_sweeps)
        .def_readonly("dense_hessian_calls", &Diagnostics::dense_hessian_calls)
        .def_readonly("branch", &Diagnostics::branch)
        .def_readonly("stable_root_count", &Diagnostics::stable_root_count);
    py::class_<epcsaft::PropertyActionBatch>(module, "PropertyActionBatch")
        .def_readonly("results", &epcsaft::PropertyActionBatch::results)
        .def_readonly("maximum_energy_order", &epcsaft::PropertyActionBatch::maximum_energy_order)
        .def_readonly("diagnostics", &epcsaft::PropertyActionBatch::diagnostics)
        .def_readonly("density_action", &epcsaft::PropertyActionBatch::density_action)
        .def_readonly("density_status", &epcsaft::PropertyActionBatch::density_status);
    py::enum_<epcsaft::ThermochemistryStatus>(module, "ThermochemistryStatus")
        .value("MissingRecord", epcsaft::ThermochemistryStatus::MissingRecord)
        .value("Available", epcsaft::ThermochemistryStatus::Available)
        .value("OutsideInterval", epcsaft::ThermochemistryStatus::OutsideInterval)
        .value("Nonfinite", epcsaft::ThermochemistryStatus::Nonfinite);
    py::class_<epcsaft::State>(module, "State")
        .def_readonly("hard_chain", &epcsaft::State::hard_chain)
        .def_readonly("dispersion", &epcsaft::State::dispersion)
        .def_readonly("association", &epcsaft::State::association)
        .def_readonly("debye_huckel", &epcsaft::State::debye_huckel)
        .def_readonly("born", &epcsaft::State::born)
        .def_readonly("bulk_relative_permittivity", &epcsaft::State::bulk_relative_permittivity)
        .def_readonly("residual_helmholtz", &epcsaft::State::residual_helmholtz)
        .def_readonly("compressibility_factor", &epcsaft::State::compressibility_factor)
        .def_readonly("pressure", &epcsaft::State::pressure)
        .def_readonly("residual_chemical_potential_over_rt", &epcsaft::State::residual_chemical_potential_over_rt)
        .def_readonly("log_fugacity_coefficient", &epcsaft::State::log_fugacity_coefficient)
        .def_readonly("residual_entropy", &epcsaft::State::residual_entropy)
        .def_readonly("residual_internal_energy", &epcsaft::State::residual_internal_energy)
        .def_readonly("residual_enthalpy", &epcsaft::State::residual_enthalpy)
        .def_readonly("residual_isochoric_heat_capacity", &epcsaft::State::residual_isochoric_heat_capacity)
        .def_readonly("residual_isobaric_heat_capacity", &epcsaft::State::residual_isobaric_heat_capacity)
        .def_readonly("ideal_enthalpy", &epcsaft::State::ideal_enthalpy)
        .def_readonly("ideal_internal_energy", &epcsaft::State::ideal_internal_energy)
        .def_readonly("ideal_entropy", &epcsaft::State::ideal_entropy)
        .def_readonly("ideal_isobaric_heat_capacity", &epcsaft::State::ideal_isobaric_heat_capacity)
        .def_readonly("ideal_isochoric_heat_capacity", &epcsaft::State::ideal_isochoric_heat_capacity)
        .def_readonly("total_enthalpy", &epcsaft::State::total_enthalpy)
        .def_readonly("total_internal_energy", &epcsaft::State::total_internal_energy)
        .def_readonly("total_entropy", &epcsaft::State::total_entropy)
        .def_readonly("total_isobaric_heat_capacity", &epcsaft::State::total_isobaric_heat_capacity)
        .def_readonly("total_isochoric_heat_capacity", &epcsaft::State::total_isochoric_heat_capacity)
        .def_readonly("molar_density", &epcsaft::State::molar_density)
        .def_readonly("pressure_density_derivative", &epcsaft::State::pressure_density_derivative)
        .def_readonly("site_fractions", &epcsaft::State::site_fractions)
        .def_readonly("density_convergence", &epcsaft::State::density_convergence)
        // Unspecified and absent in the T-rho chart: no pressure-root search was performed.
        .def_readonly("thermochemistry_status", &epcsaft::State::thermochemistry_status)
        .def_readonly("density_branch", &epcsaft::State::density_branch)
        .def_readonly("stable_root_count", &epcsaft::State::stable_root_count);
    py::class_<epcsaft::SweepCounters>(module, "SweepCounters")
        .def_readonly("forward", &epcsaft::SweepCounters::forward)
        .def_readonly("reverse", &epcsaft::SweepCounters::reverse)
        .def_readonly("hessian_calls", &epcsaft::SweepCounters::hessian_calls)
        .def_readonly("recordings", &epcsaft::SweepCounters::recordings);
    py::class_<epcsaft::HelmholtzModel>(module, "HelmholtzModel")
        .def(py::init<epcsaft::Parameters, std::vector<epcsaft::ActiveParameter>,
                      std::optional<epcsaft::ThermochemistryRecord>>(),
             py::arg("parameters"), py::arg("active_parameters") = std::vector<epcsaft::ActiveParameter>{},
             py::arg("thermochemistry") = std::nullopt)
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
        .def("state_at_pressure", &epcsaft::HelmholtzModel::state_at_pressure,
             py::arg("temperature"), py::arg("pressure"), py::arg("mole_fractions"),
             py::arg("branch"), py::arg("anchor_density") = std::nullopt)
        .def("counters", &epcsaft::HelmholtzModel::counters)
        .def("reset_counters", &epcsaft::HelmholtzModel::reset_counters);
    py::class_<epcsaft::InfiniteDilutionReference>(module, "InfiniteDilutionReference")
        .def_readonly("log_fugacity", &epcsaft::InfiniteDilutionReference::log_fugacity)
        .def_readonly("molality", &epcsaft::InfiniteDilutionReference::molality)
        .def_readonly("convergence_error", &epcsaft::InfiniteDilutionReference::convergence_error)
        .def_readonly("pure_solvent_log_fugacity", &epcsaft::InfiniteDilutionReference::pure_solvent_log_fugacity);
    py::class_<epcsaft::ElectrolyteActivity>(module, "ElectrolyteActivity")
        .def_readonly("single_ion", &epcsaft::ElectrolyteActivity::single_ion)
        .def_readonly("mean_ionic_mole_fraction_scale", &epcsaft::ElectrolyteActivity::mean_ionic_mole_fraction_scale)
        .def_readonly("mean_ionic_log_molality_scale", &epcsaft::ElectrolyteActivity::mean_ionic_log_molality_scale)
        .def_readonly("mean_ionic_molality_scale", &epcsaft::ElectrolyteActivity::mean_ionic_molality_scale)
        .def_readonly("osmotic_coefficient", &epcsaft::ElectrolyteActivity::osmotic_coefficient);
    module.def("infinite_dilution_reference", &epcsaft::infinite_dilution_reference);
    module.def("electrolyte_activity", &epcsaft::electrolyte_activity);
}
