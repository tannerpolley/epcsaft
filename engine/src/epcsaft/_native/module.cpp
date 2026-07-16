#include "eos.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(_core, module) {
    py::class_<epcsaft::NativeModel>(module, "NativeModel")
        .def(py::init<
             std::vector<double>,
             std::vector<double>,
             std::vector<double>,
             std::vector<double>>());
}
