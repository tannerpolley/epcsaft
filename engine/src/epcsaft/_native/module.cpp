#include "eos.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(_core, module) {
    py::class_<epcsaft::NativeResult>(module, "NativeResult")
        .def_readonly("hard_chain", &epcsaft::NativeResult::hard_chain)
        .def_readonly("dispersion", &epcsaft::NativeResult::dispersion)
        .def_readonly(
            "residual_helmholtz", &epcsaft::NativeResult::residual_helmholtz
        )
        .def_readonly(
            "compressibility_factor",
            &epcsaft::NativeResult::compressibility_factor
        )
        .def_readonly("pressure", &epcsaft::NativeResult::pressure);

    py::class_<epcsaft::NativeModel>(module, "NativeModel")
        .def(py::init<
             std::vector<double>,
             std::vector<double>,
             std::vector<double>,
             std::vector<double>>())
        .def("evaluate", &epcsaft::NativeModel::evaluate);
}
