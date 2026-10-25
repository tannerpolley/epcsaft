// Fixed-T,P liquid density using teqp's exact pressure derivatives.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <cmath>
#include <stdexcept>
namespace py = pybind11;

PYBIND11_MODULE(teqp_density, m) {
    m.def("pressure", [](py::object model, double T, py::array_t<double> concentrations) {
        const auto rho = concentrations.unchecked<1>();
        double total = 0.;
        for (py::ssize_t i = 0; i < rho.shape(0); ++i) total += rho(i);
        py::array_t<double> composition(rho.shape(0));
        auto x = composition.mutable_unchecked<1>();
        for (py::ssize_t i = 0; i < rho.shape(0); ++i) x(i) = rho(i)/total;
        return total*model.attr("get_R")(composition).cast<double>()*T
            + model.attr("get_pr")(T, concentrations).cast<double>();
    });
    m.def("solve", [](py::object model, double T, double p, double rho, py::object x) {
        const double R = model.attr("get_R")(x).cast<double>();
        double step = 0.;
        int iterations = 0;
        for (; iterations < 20; ++iterations) {
            const double a1 = model.attr("get_Ar01")(T, rho, x).cast<double>();
            const double a2 = model.attr("get_Ar02")(T, rho, x).cast<double>();
            const double slope = R*T*(1 + 2*a1 + a2);
            if (!(slope > 0)) throw std::runtime_error("nonpositive liquid pressure slope");
            step = (rho*R*T*(1 + a1) - p)/slope;
            rho -= step;
            if (!(rho > 0) || !std::isfinite(rho)) throw std::runtime_error("invalid liquid density");
            if (std::abs(step)/rho < 1e-14) break;
        }
        if (iterations == 20) throw std::runtime_error("density iteration limit");
        return py::make_tuple(rho, iterations + 1, std::abs(step)/rho);
    });
}
