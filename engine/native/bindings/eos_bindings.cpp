#include "model/helmholtz_model.hpp"
#include "surface/parameter_builder.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <optional>

namespace py = pybind11;

void register_eos(py::module_& module) {
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
    py::enum_<epcsaft::DensityBranch>(module, "DensityBranch")
        .value("Unspecified", epcsaft::DensityBranch::Unspecified)
        .value("Vapor", epcsaft::DensityBranch::Vapor)
        .value("Liquid", epcsaft::DensityBranch::Liquid);
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
        .def_readonly("molar_density", &epcsaft::State::molar_density)
        .def_readonly("pressure_density_derivative", &epcsaft::State::pressure_density_derivative)
        .def_readonly("site_fractions", &epcsaft::State::site_fractions)
        .def_readonly("density_convergence", &epcsaft::State::density_convergence);
    py::class_<epcsaft::SweepCounters>(module, "SweepCounters")
        .def_readonly("forward", &epcsaft::SweepCounters::forward)
        .def_readonly("reverse", &epcsaft::SweepCounters::reverse)
        .def_readonly("hessian_calls", &epcsaft::SweepCounters::hessian_calls)
        .def_readonly("recordings", &epcsaft::SweepCounters::recordings);
    py::class_<epcsaft::HelmholtzModel>(module, "HelmholtzModel")
        .def(py::init<epcsaft::Parameters>())
        .def("state", &epcsaft::HelmholtzModel::state)
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
