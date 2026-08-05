#ifndef EPCSAFT_ENGINE_REGRESSION_EVALUATOR_FIT_HPP
#define EPCSAFT_ENGINE_REGRESSION_EVALUATOR_FIT_HPP

#include "ceres_core.hpp"
#include "regression_evaluator_v1.h"

#ifndef EPCSAFT_ENGINE_REGRESSION_EVALUATOR_CORE_ONLY
#include <Python.h>
#endif

#include <cstddef>
#include <string>
#include <vector>

namespace epcsaft_engine::regression::evaluator {

struct Metadata final {
    std::string evaluator_identity;
    std::string capability_id;
    std::string capability_fingerprint;
    std::string provider_artifact_identity;
    std::string owner_artifact_identity;
    std::string contract_fingerprint;
    std::string model_fingerprint;
    std::string provider_parameter_fingerprint;
    std::string provider_topology_fingerprint;
    std::string artifact_identity;
};

struct Parameter final {
    std::string id;
    std::string unit;
    double origin;
    double scale;
    double lower;
    double upper;
};

struct Row final {
    std::string id;
    std::string partition;
    std::string state_id;
    std::string state_schema_id;
    std::string source_id;
    std::string source_locator;
    std::string primitive_id;
    std::string primitive_unit;
    std::string transform;
    std::string reference_id;
    std::string reference_fingerprint;
    double observed;
    double scale;
};

struct Problem final {
    Metadata metadata;
    std::vector<Parameter> parameters;
    std::vector<std::vector<double>> starts;
    std::vector<Row> rows;
    internal::SolverControls controls;
    double maximum_condition_number;
    double confirmation_parameter_delta;
    double confirmation_cost_delta;
};

struct RowEvaluation final {
    double value;
    double residual;
    std::vector<double> scaled_solver_jacobian;
    epcsaft_engine_regression_evaluator_row_result_v1 certificate;
};

struct FitResult final {
    std::vector<internal::SolveResult> solves;
    std::vector<RowEvaluation> rows;
    double confirmation_parameter_delta{0.0};
    double confirmation_cost_delta{0.0};
    bool confirmations_usable{false};
    std::string provider_parameter_fingerprint;
};

void validate_contract(
    const epcsaft_engine_regression_evaluator_sdk_v1& sdk,
    const Problem& problem
);

FitResult solve(
    const epcsaft_engine_regression_evaluator_sdk_v1& sdk,
    const Problem& problem
);

#ifdef EPCSAFT_ENGINE_REGRESSION_EVALUATOR_CORE_ONLY
bool evaluate_at_for_test(
    const epcsaft_engine_regression_evaluator_sdk_v1& sdk,
    const Problem& problem,
    const std::vector<double>& solver_parameters,
    std::vector<double>& residuals,
    std::vector<double>& jacobian,
    std::string& failure_reason
);
#endif

}  // namespace epcsaft_engine::regression::evaluator

namespace epcsaft_engine::regression {

#ifndef EPCSAFT_ENGINE_REGRESSION_EVALUATOR_CORE_ONLY
PyObject* solve_evaluator_python(
    PyObject* capsule, PyObject* payload
);
#endif

}  // namespace epcsaft_engine::regression

#endif
