#ifndef EPCSAFT_ENGINE_REGRESSION_PURE_SATURATION_FIT_INTERNAL_HPP
#define EPCSAFT_ENGINE_REGRESSION_PURE_SATURATION_FIT_INTERNAL_HPP

#include <Python.h>

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace epcsaft_engine::regression::internal {

constexpr std::size_t row_count = 4;
constexpr std::size_t residuals_per_row = 4;
constexpr std::size_t parameter_count = 3;
constexpr std::size_t variable_count = parameter_count + 2 * row_count;
constexpr std::size_t residual_count = residuals_per_row * row_count;

struct PyObjectDeleter final {
    void operator()(PyObject* object) const noexcept { Py_XDECREF(object); }
};
using OwnedPyObject = std::unique_ptr<PyObject, PyObjectDeleter>;

struct Row final {
    std::string row_id;
    std::string component_id;
    double temperature;
    double pressure;
    double liquid_density;
    std::string source_id;
};

struct Payload final {
    std::vector<std::string> identity;
    std::array<Row, row_count> rows;
    std::array<double, parameter_count> start;
    std::array<double, parameter_count> lower;
    std::array<double, parameter_count> upper;
    std::array<double, parameter_count> parameter_scale;
    double amount;
    double molar_mass;
    std::array<double, residuals_per_row> weights;
    std::array<double, 2> liquid_volume_bounds;
    std::array<double, 2> vapor_volume_bounds;
    double topology_separation;
    int max_iterations;
    double function_tolerance;
    double gradient_tolerance;
    double parameter_tolerance;
    std::array<double, 2> reporting_pressure_bounds;
    double confirmation_liquid_start_multiplier;
    double confirmation_vapor_start_multiplier;
    double confirmation_parameter_delta;
    double confirmation_cost_delta;
    double reporting_pressure_closure;
    double reporting_mu_closure;
    int num_threads;
};

bool positive_finite(double value);
std::vector<double> doubles(PyObject* object, std::size_t expected, const char* label);
Payload parse_payload(PyObject* object);
Row parse_row(PyObject* object, const std::string& component_id, std::size_t source_index);
std::size_t reporting_row_count(const Payload& payload);

}  // namespace epcsaft_engine::regression::internal

#endif
