#include "eos.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <utility>
#include <vector>

namespace py = pybind11;

PYBIND11_MODULE(_core, module) {
    py::class_<epcsaft::NativeResult>(module, "NativeResult")
        .def_readonly("hard_chain", &epcsaft::NativeResult::hard_chain)
        .def_readonly("dispersion", &epcsaft::NativeResult::dispersion)
        .def_readonly("association", &epcsaft::NativeResult::association)
        .def_readonly("debye_huckel", &epcsaft::NativeResult::debye_huckel)
        .def_readonly(
            "residual_helmholtz", &epcsaft::NativeResult::residual_helmholtz
        )
        .def_readonly(
            "compressibility_factor",
            &epcsaft::NativeResult::compressibility_factor
        )
        .def_readonly("pressure", &epcsaft::NativeResult::pressure);

    py::class_<epcsaft::NativeModel>(module, "NativeModel")
        .def(py::init([](
            std::vector<double> segment_count,
            std::vector<double> sigma_constant,
            std::vector<int> sigma_term_component,
            std::vector<double> sigma_term_amplitude,
            std::vector<double> sigma_term_exponent,
            std::vector<double> dispersion_energy_over_k,
            std::vector<double> k_ij,
            std::vector<int> site_component,
            std::vector<int> site_multiplicity,
            std::vector<int> association_active,
            std::vector<double> association_energy_over_k,
            std::vector<double> association_volume,
            std::vector<double> molar_mass_kg_per_mol,
            std::vector<int> charge_number,
            std::vector<double> relative_permittivity,
            double dielectric_ion_suppression_coefficient
        ) {
            return epcsaft::NativeModel(epcsaft::NativeParameters{
                std::move(segment_count),
                std::move(sigma_constant),
                std::move(sigma_term_component),
                std::move(sigma_term_amplitude),
                std::move(sigma_term_exponent),
                std::move(dispersion_energy_over_k),
                std::move(k_ij),
                std::move(site_component),
                std::move(site_multiplicity),
                std::move(association_active),
                std::move(association_energy_over_k),
                std::move(association_volume),
                std::move(molar_mass_kg_per_mol),
                std::move(charge_number),
                std::move(relative_permittivity),
                dielectric_ion_suppression_coefficient,
            });
        }))
        .def("evaluate", &epcsaft::NativeModel::evaluate);
}
