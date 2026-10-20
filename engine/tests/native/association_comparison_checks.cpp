#include <cppad/example/cppad_eigen.hpp>

#include "contributions/association.hpp"
#include "contributions/helmholtz.hpp"
#include "properties/helmholtz_tape.hpp"
#include "properties_support.hpp"
#include "recording_support.hpp"

#include <Eigen/LU>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#define CHECK(condition) recording_support::check((condition), __FILE__, __LINE__, #condition)

namespace {

using AD = CppAD::AD<double>;

// The association-comparison analysis reads these expected keys and group counts
// (EPCSAFT_EXPECT_JSON) to confirm that every comparison it needs was recorded.
void expect_keys() {
    std::printf("EPCSAFT_EXPECT_JSON {\"schema_version\":1,\"test_id\":\"native.association.comparison\","
                "\"expected_comparison_keys\":[\"site-fraction/*\",\"candidate-residual/*\",\"reference-residual/*\","
                "\"association-energy/*\",\"residual-energy/*\",\"association-gradient/*\",\"residual-gradient/*\","
                "\"association-hessian/*\",\"residual-hessian/*\",\"reference-quality/*\",\"centered-gradient/*\","
                "\"centered-hessian/*\",\"centered-mixed/*\",\"analytic-control/*\",\"k-residual-identity/*\","
                "\"production-weighted-hessian/*\",\"extensive-scaling/*\",\"failure-continuation\"]}\n");
}

void expect_group(const char* test_id, const char* comparison_group_id, const char* comparison_id_prefix,
                  std::size_t expected_count, const std::string& definition, const char* coordinates) {
    using recording_support::json_string;
    std::printf("EPCSAFT_EXPECT_JSON {\"schema_version\":1,\"test_id\":%s,\"comparison_group_id\":%s,"
                "\"comparison_id_prefix\":%s,\"expected_count\":%zu,\"coordinates\":%s,\"definition\":%s}\n",
                json_string(test_id).c_str(), json_string(comparison_group_id).c_str(),
                json_string(comparison_id_prefix).c_str(), expected_count,
                recording_support::json_fields_object(coordinates).c_str(),
                recording_support::json_fields_object(definition).c_str());
    std::fflush(stdout);
}

using Parameters = epcsaft::Parameters;
constexpr std::array<std::size_t, 4> DEPTHS = {8, 16, 24, 32};
constexpr std::size_t QUALIFYING_DEPTH = 16;
constexpr double VALUE_TOL = 1.0e-10;
constexpr double SITE_RTOL = 1.0e-8;
constexpr double ENERGY_RTOL = 1.0e-9;
constexpr double REFERENCE_TOL = 1.0e-12;
constexpr double CONDITION_LIMIT = 1.0e6;
constexpr double DERIVATIVE_ATOL = 1.0e-8;
constexpr double DERIVATIVE_RTOL = 1.0e-6;
constexpr double CENTERED_ATOL = 1.0e-7;
constexpr double CENTERED_RTOL = 1.0e-5;
constexpr double SENSITIVITY_RESIDUAL_LIMIT = 1.0e-10;
constexpr double SENSITIVITY_PROXY_LIMIT = 1.0e-8;
constexpr double NEWTON_DEFECT_LIMIT = 1.0e-10;

// #705: captured rejected trial at 298.15 K; density in mol/m^3, strengths in Angstrom^3.
// The previous step bound exhausted 50 iterations at residual 5.480282894154698e-12.
void association_unit_ceiling() {
    constexpr double temperature = 298.15;
    constexpr double density = 52992.32016447701;
    const std::vector<double> mole_fractions = {
        0.9967664629120051, 9.581789719163911e-29, 0.003232311646483076,
        6.310494279431576e-07, 5.712187821800883e-10, 5.938208650913641e-07};
    Parameters parameters;
    parameters.site_component = {0, 0, 1, 1, 2, 2, 3, 3};
    parameters.site_multiplicity = std::vector<int>(8, 1);
    const std::vector<double> strength = {
        0.0, 15964.491628191587, 0.0, 1986249.8013975348, 0.0, 4550137.1076187035, 0.0, 0.0,
        15964.491628191587, 0.0, 1986249.8013975348, 0.0, 4550137.1076187035, 0.0, 0.0, 0.0,
        0.0, 1986249.8013975348, 0.0, 286578809.8455787, 0.0, 633297185.5581268, 0.0, 61842.874171620366,
        1986249.8013975348, 0.0, 286578809.8455787, 0.0, 633297185.5581268, 0.0, 61842.874171620366, 0.0,
        0.0, 4550137.1076187035, 0.0, 633297185.5581268, 0.0, 1412151858.732317, 0.0, 0.0,
        4550137.1076187035, 0.0, 633297185.5581268, 0.0, 1412151858.732317, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 61842.874171620366, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 61842.874171620366, 0.0, 0.0, 0.0, 0.0, 0.0,
    };
    epcsaft::core::Convergence convergence;
    const auto fractions = epcsaft::solve_site_fractions(temperature, density, mole_fractions,
        parameters, {}, {}, &convergence, nullptr, &strength);
    const auto residuals = epcsaft::association_residuals(density, mole_fractions,
        parameters, fractions, strength);
    double residual = 0.0;
    for (double value : residuals) {
        CHECK(std::isfinite(value));
        residual = std::max(residual, std::abs(value));
    }
    const double tolerance = epcsaft::policy::SolverOptions{}.value(epcsaft::policy::association);
    CHECK(residual <= tolerance);
    CHECK(convergence.iterations <= epcsaft::policy::ASSOCIATION_MAX_ITERATIONS.value);
    CHECK(std::all_of(fractions.begin(), fractions.end(),
        [](double value) { return std::isfinite(value) && value > 0.0 && value <= 1.0; }));
    std::printf("association unit ceiling: residual=%.17g tolerance=%.17g iterations=%zu\n",
        residual, tolerance, convergence.iterations);
}

void cppad_error_handler(bool, int, const char*, const char*, const char* message) {
    throw std::runtime_error(std::string("CppAD error: ") + message);
}
struct TestCase final {
    const char* id;
    Parameters parameters;
    double temperature;
    double volume;
    std::vector<double> amounts;
    std::optional<epcsaft::ActiveParameter> active = std::nullopt;
    bool analytic = false;
    std::optional<std::uint32_t> seed = std::nullopt;
    const char* partition = "controls";
    const char* stratum = "controls";
};
struct Tape final {
    CppAD::ADFun<double> function;
    std::vector<double> point;
    std::size_t width = 0;
    std::size_t sites = 0;
};
struct Candidate final {
    double association = std::numeric_limits<double>::quiet_NaN();
    double residual = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> sites;
    std::vector<double> mass_action;
    std::vector<double> association_gradient;
    std::vector<double> association_hessian;
    std::vector<double> residual_gradient;
    std::vector<double> residual_hessian;
};
struct EnergyValues final {
    double association = std::numeric_limits<double>::quiet_NaN();
    double residual = std::numeric_limits<double>::quiet_NaN();
    bool available = false;
};
struct Reference final {
    std::vector<double> sites;
    std::vector<double> mass_action;
    double association = std::numeric_limits<double>::quiet_NaN();
    double residual = std::numeric_limits<double>::quiet_NaN();
    properties_support::ReferenceTensor tensor;
    properties_support::ReferenceTensor association_tensor;
    double condition = std::numeric_limits<double>::quiet_NaN();
    std::string status = "unavailable";
    std::string value_status = "unavailable";
    std::string derivative_status = "unavailable";
    std::string error;
    double newton_defect = std::numeric_limits<double>::quiet_NaN();
    bool value_available = false;
    bool derivative_available = false;
};

std::vector<epcsaft::ActiveParameter> active_parameters(const TestCase& test_case) {
    return test_case.active.has_value()
        ? std::vector<epcsaft::ActiveParameter>{*test_case.active}
        : std::vector<epcsaft::ActiveParameter>{};
}

std::vector<double> composition(const std::vector<double>& amounts) {
    double total = 0.0;
    for (double amount : amounts) total += amount;
    std::vector<double> result;
    result.reserve(amounts.size());
    for (double amount : amounts) result.push_back(amount / total);
    return result;
}
template <typename Scalar, typename ParameterSet>
std::vector<Scalar> association_k(const Scalar& density, const std::vector<Scalar>& composition,
    const ParameterSet& parameters, const std::vector<Scalar>& delta) {
    const std::size_t sites = parameters.site_component.size();
    const Scalar number_density = density * epcsaft::NUMBER_DENSITY_PER_ANGSTROM3;
    std::vector<Scalar> result(sites * sites, Scalar(0.0));
    for (std::size_t a = 0; a < sites; ++a) for (std::size_t b = 0; b < sites; ++b) {
        const auto component = static_cast<std::size_t>(parameters.site_component[b]);
        result[a * sites + b] = number_density * composition[component]
            * parameters.site_multiplicity[b] * delta[a * sites + b];
    }
    return result;
}
template <typename Scalar>
std::vector<Scalar> solve_unpivoted(std::vector<Scalar> matrix, std::vector<Scalar> rhs) {
    const std::size_t size = rhs.size();
    for (std::size_t pivot = 0; pivot < size; ++pivot) {
        const Scalar diagonal = matrix[pivot * size + pivot];
        for (std::size_t row = pivot + 1; row < size; ++row) {
            const Scalar factor = matrix[row * size + pivot] / diagonal;
            for (std::size_t column = pivot + 1; column < size; ++column)
                matrix[row * size + column] -= factor * matrix[pivot * size + column];
            rhs[row] -= factor * rhs[pivot];
        }
    }
    std::vector<Scalar> result(size, Scalar(0.0));
    for (std::size_t row = size; row-- > 0;) {
        Scalar value = rhs[row];
        for (std::size_t column = row + 1; column < size; ++column)
            value -= matrix[row * size + column] * result[column];
        result[row] = value / matrix[row * size + row];
    }
    return result;
}
template <typename Scalar, typename ParameterSet>
std::vector<Scalar> fixed_candidate(const Scalar& density, const std::vector<Scalar>& composition,
    const ParameterSet& parameters, const std::vector<Scalar>& delta, std::size_t depth) {
    const std::size_t sites = parameters.site_component.size();
    const auto k = association_k(density, composition, parameters, delta);
    std::vector<Scalar> y(sites, Scalar(0.0));
    for (std::size_t row = 0; row < sites; ++row) {
        Scalar row_sum = Scalar(0.0);
        for (std::size_t column = 0; column < sites; ++column) row_sum += k[row * sites + column];
        y[row] = epcsaft::scalar_log(Scalar(2.0) / (Scalar(1.0)
            + epcsaft::scalar_sqrt(Scalar(1.0) + Scalar(4.0) * row_sum)));
    }
    for (std::size_t iteration = 0; iteration < depth; ++iteration) {
        std::vector<Scalar> x(sites, Scalar(0.0));
        for (std::size_t row = 0; row < sites; ++row) x[row] = epcsaft::scalar_exp(y[row]);
        std::vector<Scalar> bonded(sites, Scalar(0.0));
        for (std::size_t row = 0; row < sites; ++row)
            for (std::size_t column = 0; column < sites; ++column)
                bonded[row] += k[row * sites + column] * x[column];
        std::vector<Scalar> residual(sites, Scalar(0.0));
        std::vector<Scalar> jacobian(sites * sites, Scalar(0.0));
        for (std::size_t row = 0; row < sites; ++row) {
            const Scalar denominator = Scalar(1.0) + bonded[row];
            residual[row] = y[row] + epcsaft::scalar_log(denominator);
            for (std::size_t column = 0; column < sites; ++column)
                jacobian[row * sites + column] = (row == column ? Scalar(1.0) : Scalar(0.0))
                    + k[row * sites + column] * x[column] / denominator;
        }
        const auto step = solve_unpivoted(jacobian, residual);
        Scalar squared_norm = Scalar(0.0);
        for (const Scalar& value : step) squared_norm += value * value;
        const Scalar limiter = epcsaft::scalar_sqrt(Scalar(1.0) + squared_norm / Scalar(64.0));
        for (std::size_t row = 0; row < sites; ++row) y[row] -= step[row] / limiter;
    }
    std::vector<Scalar> result(sites, Scalar(0.0));
    for (std::size_t site = 0; site < sites; ++site) result[site] = epcsaft::scalar_exp(y[site]);
    return result;
}
std::vector<double> point_for(const TestCase& test_case) {
    std::vector<double> point = {test_case.temperature, test_case.volume};
    point.insert(point.end(), test_case.amounts.begin(), test_case.amounts.end());
    if (test_case.active.has_value()) {
        const auto& active = *test_case.active;
        const auto& values = active.family == epcsaft::ParameterFamily::AssociationEnergyOverK
            ? test_case.parameters.association_energy_over_k : test_case.parameters.association_volume;
        point.push_back(values[active.indices.front()]);
    }
    return point;
}
Tape make_tape(const TestCase& test_case, std::size_t depth) {
    const std::size_t components = test_case.parameters.component_count();
    const std::vector<double> point = point_for(test_case);
    std::vector<AD> independent(point.begin(), point.end());
    CppAD::Independent(independent);
    std::vector<AD> amounts(independent.begin() + 2,
        independent.begin() + 2 + static_cast<std::ptrdiff_t>(components));
    AD total = AD(0.0);
    for (const AD& amount : amounts) total += amount;
    const AD density = total / independent[1];
    std::vector<AD> mole_fractions;
    for (const AD& amount : amounts) mole_fractions.push_back(amount / total);
    std::vector<epcsaft::detail::ActiveParameter> active;
    if (test_case.active.has_value()) active.push_back(*test_case.active);
    const auto parameters = epcsaft::detail::make_tape_parameters(test_case.parameters, active,
        std::vector<AD>(independent.begin() + 2 + static_cast<std::ptrdiff_t>(components), independent.end()));
    const auto sigma = epcsaft::segment_sigma(independent[0], parameters);
    const auto diameter = epcsaft::segment_diameter(independent[0], parameters, sigma);
    const auto delta = epcsaft::association_delta_values(independent[0], density, mole_fractions,
        parameters, sigma, diameter, parameters.association_energy_over_k, parameters.association_volume);
    const auto sites = fixed_candidate(density, mole_fractions, parameters, delta, depth);
    const auto terms = epcsaft::residual_helmholtz(independent[0], density, mole_fractions,
        parameters, sigma, diameter, sites, AD(0.0), AD(0.0));
    const AD association = total * terms.association;
    const AD residual = total * (terms.hard_chain + terms.dispersion + terms.association);
    const auto mass_action = epcsaft::association_residuals(density, mole_fractions, parameters, sites, delta);
    std::vector<AD> dependent = {association, residual};
    dependent.insert(dependent.end(), sites.begin(), sites.end());
    dependent.insert(dependent.end(), mass_action.begin(), mass_action.end());
    Tape result;
    result.function = CppAD::ADFun<double>(independent, dependent);
    result.point = point;
    result.width = point.size();
    result.sites = sites.size();
    return result;
}
Candidate evaluate_candidate(const TestCase& test_case, std::size_t depth) {
    Tape tape = make_tape(test_case, depth);
    const auto values = tape.function.Forward(0, tape.point);
    const auto jacobian = tape.function.Jacobian(tape.point);
    Candidate result;
    result.association = values[0];
    result.residual = values[1];
    result.sites.assign(values.begin() + 2, values.begin() + 2 + tape.sites);
    result.mass_action.assign(values.begin() + 2 + tape.sites, values.end());
    result.association_gradient.resize(tape.width);
    result.residual_gradient.resize(tape.width);
    for (std::size_t coordinate = 0; coordinate < tape.width; ++coordinate) {
        result.association_gradient[coordinate] = jacobian[coordinate];
        result.residual_gradient[coordinate] = jacobian[tape.width + coordinate];
    }
    result.association_hessian = tape.function.Hessian(tape.point, 0);
    result.residual_hessian = tape.function.Hessian(tape.point, 1);
    return result;
}
Reference evaluate_reference(const TestCase& test_case, bool derivatives = true) {
    Reference result;
    double total = 0.0;
    for (double amount : test_case.amounts) total += amount;
    const double density = total / test_case.volume;
    const auto mole_fractions = composition(test_case.amounts);
    const auto sigma = epcsaft::segment_sigma(test_case.temperature, test_case.parameters);
    const auto diameter = epcsaft::segment_diameter(test_case.temperature, test_case.parameters, sigma);
    result.sites = epcsaft::solve_site_fractions(test_case.temperature, density, mole_fractions,
        test_case.parameters, sigma, diameter);
    const auto delta = epcsaft::association_delta(test_case.temperature, density, mole_fractions,
        test_case.parameters, sigma, diameter);
    result.mass_action = epcsaft::association_residuals(density, mole_fractions, test_case.parameters,
        result.sites, delta);
    const auto terms = epcsaft::residual_helmholtz(test_case.temperature, density, mole_fractions,
        test_case.parameters, sigma, diameter, result.sites, 0.0, 0.0);
    result.association = total * terms.association;
    result.residual = total * (terms.hard_chain + terms.dispersion + terms.association);
    const std::size_t sites = result.sites.size();
    const bool finite_sites = std::all_of(result.sites.begin(), result.sites.end(),
        [](double value) { return std::isfinite(value) && value > 0.0 && value <= 1.0; });
    double residual_norm = 0.0;
    bool finite_residual = true;
    for (double value : result.mass_action) {
        finite_residual = finite_residual && std::isfinite(value);
        if (std::isfinite(value)) residual_norm = std::max(residual_norm, std::abs(value));
    }
    bool value_quality = finite_sites && finite_residual && std::isfinite(result.association)
        && std::isfinite(result.residual) && residual_norm <= REFERENCE_TOL;
    if (!value_quality) result.value_status = "reference_value_failure";
    if (value_quality && sites != 0) {
        const auto k = association_k(density, mole_fractions, test_case.parameters, delta);
        Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(sites), static_cast<Eigen::Index>(sites));
        for (std::size_t row = 0; row < sites; ++row) {
            double bonded = 0.0;
            for (std::size_t column = 0; column < sites; ++column)
                bonded += k[row * sites + column] * result.sites[column];
            for (std::size_t column = 0; column < sites; ++column)
                jacobian(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) =
                    (row == column ? 1.0 + bonded : 0.0) + result.sites[row] * k[row * sites + column];
        }
        const Eigen::FullPivLU<Eigen::MatrixXd> factorization(jacobian);
        if (factorization.rank() != jacobian.rows()) {
            value_quality = false;
            result.value_status = "reference_rank_failure";
        } else {
            const Eigen::VectorXd residual = Eigen::Map<const Eigen::VectorXd>(
                result.mass_action.data(), static_cast<Eigen::Index>(sites));
            const Eigen::VectorXd step = factorization.solve(residual);
            result.newton_defect = 0.0;
            for (std::size_t site = 0; site < sites; ++site) {
                const double defect = step(static_cast<Eigen::Index>(site)) / result.sites[site];
                if (!std::isfinite(defect)) result.newton_defect = std::numeric_limits<double>::quiet_NaN();
                else result.newton_defect = std::max(result.newton_defect, std::abs(defect));
            }
            value_quality = std::isfinite(result.newton_defect)
                && result.newton_defect <= NEWTON_DEFECT_LIMIT;
            if (!value_quality) result.value_status = "reference_newton_defect";
        }
    } else if (value_quality) {
        result.newton_defect = 0.0;
    }
    result.value_available = value_quality;
    if (value_quality && derivatives) {
        try {
            result.tensor = properties_support::reference_tensor(test_case.parameters, test_case.temperature,
                test_case.volume, test_case.amounts, result.sites, active_parameters(test_case));
            result.association_tensor = properties_support::reference_tensor(test_case.parameters,
                test_case.temperature, test_case.volume, test_case.amounts, result.sites,
                active_parameters(test_case), true);
            result.condition = result.tensor.scaled_condition;
            result.value_status = result.value_available ? "passed" : "reference_newton_defect";
            const auto finite_tensor = [](const properties_support::ReferenceTensor& tensor) {
                return std::all_of(tensor.gradient.begin(), tensor.gradient.end(),
                    [](double value) { return std::isfinite(value); })
                    && std::all_of(tensor.hessian.begin(), tensor.hessian.end(),
                        [](double value) { return std::isfinite(value); })
                    && std::isfinite(tensor.scaled_condition)
                    && std::isfinite(tensor.first_sensitivity_residual)
                    && std::isfinite(tensor.second_sensitivity_residual)
                    && std::isfinite(tensor.sensitivity_error_proxy);
            };
            const auto qualified_tensor = [&](const properties_support::ReferenceTensor& tensor) {
                return finite_tensor(tensor) && tensor.scaled_condition <= CONDITION_LIMIT
                    && tensor.first_sensitivity_residual <= SENSITIVITY_RESIDUAL_LIMIT
                    && tensor.second_sensitivity_residual <= SENSITIVITY_RESIDUAL_LIMIT
                    && tensor.sensitivity_error_proxy <= SENSITIVITY_PROXY_LIMIT;
            };
            result.derivative_available = qualified_tensor(result.tensor)
                && qualified_tensor(result.association_tensor);
            result.derivative_status = result.derivative_available ? "passed" : "reference_derivative_uncertain";
            if (result.value_available && result.derivative_available) result.status = "passed";
            else if (result.value_available) result.status = result.derivative_status;
            else result.status = result.value_status;
        } catch (const std::exception& error) {
            AD::abort_recording();
            result.error = error.what();
            result.value_status = result.value_available ? "passed" : "reference_tensor_failure";
            result.derivative_status = "reference_tensor_failure";
            result.status = result.value_available ? result.derivative_status : result.value_status;
        }
    }
    if (!result.value_available && result.value_status == "unavailable")
        result.value_status = "reference_value_failure";
    if (!result.derivative_available && result.derivative_status == "unavailable")
        result.derivative_status = "reference_derivative_unavailable";
    return result;
}
Parameters pure_two_b() {
    auto result = properties_support::ethane();
    result.site_component = {0, 0};
    result.site_multiplicity = {1, 1};
    result.association_active = {0, 1, 1, 0};
    result.association_energy_over_k = {0.0, 1800.0, 1800.0, 0.0};
    result.association_volume = {0.0, 0.10, 0.10, 0.0};
    return result;
}
Parameters diluted_two_b() {
    auto result = properties_support::neutral_binary();
    result.site_component = {0, 0};
    result.site_multiplicity = {1, 1};
    result.association_active = {0, 1, 1, 0};
    result.association_energy_over_k = {0.0, 1800.0, 1800.0, 0.0};
    result.association_volume = {0.0, 0.10, 0.10, 0.0};
    return result;
}
Parameters pure_three_b() {
    auto result = properties_support::ethane();
    result.site_component = {0, 0, 0};
    result.site_multiplicity = {1, 1, 1};
    result.association_active.assign(9, 0);
    result.association_energy_over_k.assign(9, 0.0);
    result.association_volume.assign(9, 0.0);
    for (const auto edge : {std::pair<std::size_t, std::size_t>{0, 1}, {0, 2}}) {
        result.association_active[edge.first * 3 + edge.second] = result.association_active[edge.second * 3 + edge.first] = 1;
        result.association_energy_over_k[edge.first * 3 + edge.second] = result.association_energy_over_k[edge.second * 3 + edge.first] = 1500.0;
        result.association_volume[edge.first * 3 + edge.second] = result.association_volume[edge.second * 3 + edge.first] = 0.10;
    }
    return result;
}
Parameters disconnected() {
    auto result = properties_support::neutral_binary();
    result.site_component = {0, 0, 1, 1};
    result.site_multiplicity = {1, 1, 1, 1};
    result.association_active.assign(16, 0);
    result.association_energy_over_k.assign(16, 0.0);
    result.association_volume.assign(16, 0.0);
    result.association_active[1] = result.association_active[4] = 1;
    result.association_energy_over_k[1] = result.association_energy_over_k[4] = 2200.0;
    result.association_volume[1] = result.association_volume[4] = 0.12;
    return result;
}
Parameters neutral_three_component() {
    auto result = properties_support::neutral_ternary();
    result.site_component = {0, 0, 1, 1, 2, 2};
    result.site_multiplicity = {2, 1, 1, 2, 1, 1};
    result.association_active.assign(36, 0);
    result.association_energy_over_k.assign(36, 0.0);
    result.association_volume.assign(36, 0.0);
    for (const auto edge : {std::pair<std::size_t, std::size_t>{0, 1}, {0, 3}, {2, 5}, {1, 4}}) {
        result.association_active[edge.first * 6 + edge.second] = result.association_active[edge.second * 6 + edge.first] = 1;
        result.association_energy_over_k[edge.first * 6 + edge.second] = result.association_energy_over_k[edge.second * 6 + edge.first] = 2200.0;
        result.association_volume[edge.first * 6 + edge.second] = result.association_volume[edge.second * 6 + edge.first] = 0.12;
    }
    return result;
}

TestCase tuned(const char* id, Parameters parameters, double temperature, double volume,
    std::vector<double> amounts, std::optional<epcsaft::ActiveParameter> active = std::nullopt,
    bool analytic = false, std::optional<std::uint32_t> seed = std::nullopt,
    const char* partition = "controls", const char* stratum = "controls") {
    return {id, std::move(parameters), temperature, volume, std::move(amounts), active, analytic,
        seed, partition, stratum};
}

TestCase seeded_case(const char* id, std::uint32_t seed) {
    auto result = properties_support::neutral_ternary();
    result.site_component = {0, 0, 1, 1, 2, 2};
    result.site_multiplicity = {2, 1, 1, 2, 1, 1};
    result.association_active.assign(36, 0);
    result.association_energy_over_k.assign(36, 0.0);
    result.association_volume.assign(36, 0.0);
    const std::vector<std::pair<std::size_t, std::size_t>> edges = seed % 2 == 0
        ? std::vector<std::pair<std::size_t, std::size_t>>{{0, 1}, {0, 3}, {2, 5}, {1, 4}}
        : std::vector<std::pair<std::size_t, std::size_t>>{{0, 2}, {0, 4}, {2, 3}, {1, 5}};
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> energy(1500.0, 2400.0);
    std::uniform_real_distribution<double> volume(0.08, 0.14);
    for (const auto edge : edges) {
        const double epsilon = energy(generator);
        const double delta = volume(generator);
        result.association_active[edge.first * 6 + edge.second] = result.association_active[edge.second * 6 + edge.first] = 1;
        result.association_energy_over_k[edge.first * 6 + edge.second] = result.association_energy_over_k[edge.second * 6 + edge.first] = epsilon;
        result.association_volume[edge.first * 6 + edge.second] = result.association_volume[edge.second * 6 + edge.first] = delta;
    }
    std::uniform_real_distribution<double> amount(0.1, 0.9);
    return tuned(id, std::move(result), 330.0 + 80.0 * amount(generator), 7.0e-5 + 8.0e-5 * amount(generator),
        {amount(generator), amount(generator), amount(generator)}, std::nullopt, false, seed);
}

struct GeneratedStratum final {
    const char* name;
    std::size_t components;
    std::size_t sites;
    double eta_low;
    double eta_high;
    int topology;
};

constexpr std::array<GeneratedStratum, 8> GENERATED_STRATA = {{
    {"ordinary-dilute", 1, 2, 1.0e-8, 1.0e-4, 0},
    {"ordinary-dense", 2, 6, 0.05, 0.20, 1},
    {"strong-dense", 3, 24, 0.20, 0.50, 3},
    {"heterogeneous-cross-capacity", 5, 12, 0.01, 0.15, 2},
    {"disconnected-nonassociating", 8, 16, 1.0e-4, 0.08, 0},
    {"trace-support", 2, 8, 1.0e-6, 0.10, 2},
    {"exact-zero-support", 3, 6, 0.01, 0.20, 1},
    {"parameter-boundaries", 5, 20, 1.0e-8, 0.50, 3},
}};

double generated_unit(std::mt19937& generator) {
    return static_cast<double>(generator()) / static_cast<double>(generator.max());
}

double generated_between(std::mt19937& generator, double low, double high) {
    return low + (high - low) * generated_unit(generator);
}

Parameters generated_parameters(const GeneratedStratum& stratum, std::size_t index,
    std::mt19937& generator) {
    auto result = properties_support::parameters(stratum.components);
    for (std::size_t component = 0; component < stratum.components; ++component) {
        result.segment_count[component] = 1.0 + static_cast<double>((component + index) % 5);
        const std::size_t endpoint_slot = index % 12;
        result.sigma_constant[component] = component == 0 && endpoint_slot == 0 ? 2.0
            : component == 0 && endpoint_slot == 1 ? 4.0
            : generated_between(generator, 2.0, 4.0);
        result.dispersion_energy_over_k[component] = component == 0 && endpoint_slot == 0 ? 100.0
            : component == 0 && endpoint_slot == 1 ? 350.0
            : generated_between(generator, 100.0, 350.0);
        result.molar_mass_kg_per_mol[component] = generated_between(generator, 0.015, 0.18);
    }
    result.site_component.resize(stratum.sites);
    result.site_multiplicity.resize(stratum.sites);
    for (std::size_t site = 0; site < stratum.sites; ++site) {
        result.site_component[site] = static_cast<int>(site % stratum.components);
        result.site_multiplicity[site] = 1 + static_cast<int>((site + index) % 4);
    }
    result.association_active.assign(stratum.sites * stratum.sites, 0);
    result.association_energy_over_k.assign(stratum.sites * stratum.sites, 0.0);
    result.association_volume.assign(stratum.sites * stratum.sites, 0.0);
    for (std::size_t row = 0; row < stratum.sites; ++row) for (std::size_t column = 0; column < stratum.sites; ++column) {
        if (row == column) continue;
        const std::size_t group_row = row / (stratum.topology == 1 ? 3 : 2);
        const std::size_t group_column = column / (stratum.topology == 1 ? 3 : 2);
        bool active = false;
        if (stratum.topology == 0) active = group_row == group_column && row % 2 != column % 2;
        if (stratum.topology == 1) active = group_row == group_column && (row % 3 == 0 || column % 3 == 0);
        if (stratum.topology == 2) active = (row + column + index) % 5 != 0;
        if (stratum.topology == 3) active = (row + column + index) % 3 != 0;
        if (!active) continue;
        const std::size_t pair = row * stratum.sites + column;
        result.association_active[pair] = 1;
        result.association_energy_over_k[pair] = generated_between(generator,
            stratum.topology == 3 ? 2500.0 : 100.0, stratum.topology == 3 ? 6000.0 : 3500.0);
        result.association_volume[pair] = generated_between(generator, 0.005, 0.20);
        if (stratum.name == std::string("parameter-boundaries") && (row + column + index) % 11 == 0) {
            result.association_energy_over_k[pair] = (index % 2 == 0) ? 0.0 : 6000.0;
            result.association_volume[pair] = (index % 3 == 0) ? 0.0 : 0.20;
        }
    }
    for (std::size_t row = 0; row < stratum.sites; ++row) for (std::size_t column = row + 1; column < stratum.sites; ++column) {
        const std::size_t left = row * stratum.sites + column;
        const std::size_t right = column * stratum.sites + row;
        result.association_active[right] = result.association_active[left];
        result.association_energy_over_k[right] = result.association_energy_over_k[left];
        result.association_volume[right] = result.association_volume[left];
    }
    if (stratum.topology == 0 && stratum.components > 1) {
        const int nonassociating_component = static_cast<int>(stratum.components - 1);
        for (std::size_t row = 0; row < stratum.sites; ++row) for (std::size_t column = 0; column < stratum.sites; ++column)
            if (result.site_component[row] == nonassociating_component
                || result.site_component[column] == nonassociating_component) {
                const std::size_t pair = row * stratum.sites + column;
                result.association_active[pair] = 0;
                result.association_energy_over_k[pair] = 0.0;
                result.association_volume[pair] = 0.0;
            }
    }
    return result;
}

std::array<std::string, 96> generated_names;

TestCase generated_case(std::size_t stratum_index, std::size_t repeat, bool heldout, std::size_t name_index) {
    const auto& stratum = GENERATED_STRATA[stratum_index];
    const std::uint32_t seed = 21703000u + static_cast<std::uint32_t>(name_index * 37 + (heldout ? 100000 : 0));
    std::mt19937 generator(seed);
    auto parameters = generated_parameters(stratum, name_index, generator);
    const double temperature = repeat == 0 && stratum_index == 7 ? 200.0
        : repeat == 1 && stratum_index == 7 ? 650.0 : generated_between(generator, 200.0, 650.0);
    const double total = generated_between(generator, 0.1, 2.0);
    std::vector<double> amounts(stratum.components, total / static_cast<double>(stratum.components));
    for (double& amount : amounts) amount *= generated_between(generator, 0.5, 1.5);
    if (stratum_index == 5) amounts[0] = total * 1.0e-12;
    if (stratum_index == 6) amounts[0] = 0.0;
    double amount_sum = 0.0;
    for (double amount : amounts) amount_sum += amount;
    if (amount_sum == 0.0) amounts[0] = total;
    const auto fractions = composition(amounts);
    const auto sigma = epcsaft::segment_sigma(temperature, parameters);
    const auto diameter = epcsaft::segment_diameter(temperature, parameters, sigma);
    const double eta = repeat == 0 && stratum_index == 0 ? 1.0e-8
        : repeat == 1 && stratum_index == 7 ? 0.50
        : generated_between(generator, stratum.eta_low, stratum.eta_high);
    double moment = 0.0;
    for (std::size_t component = 0; component < stratum.components; ++component)
        moment += fractions[component] * parameters.segment_count[component]
            * diameter[component] * diameter[component] * diameter[component];
    const double density = eta / (epcsaft::PI / 6.0 * epcsaft::NUMBER_DENSITY_PER_ANGSTROM3 * moment);
    const double volume = amount_sum / density;
    generated_names[name_index] = std::string(heldout ? "heldout-" : "development-")
        + stratum.name + "-" + std::to_string(repeat);
    return tuned(generated_names[name_index].c_str(), std::move(parameters), temperature, volume,
        std::move(amounts), std::nullopt, false, seed, heldout ? "heldout" : "development", stratum.name);
}

std::vector<TestCase> generated_cases() {
    std::vector<TestCase> result;
    result.reserve(96);
    std::size_t name_index = 0;
    for (std::size_t stratum = 0; stratum < GENERATED_STRATA.size(); ++stratum) {
        for (std::size_t repeat = 0; repeat < 8; ++repeat)
            result.push_back(generated_case(stratum, repeat, false, name_index++));
        for (std::size_t repeat = 0; repeat < 4; ++repeat)
            result.push_back(generated_case(stratum, repeat, true, name_index++));
    }
    return result;
}

std::string campaign_selector() {
    const char* value = std::getenv("EPCSAFT_ASSOCIATION_SELECTOR");
    return value == nullptr ? "all" : value;
}

bool selected_case(const TestCase& test_case, const std::string& selector) {
    if (selector == "all") return true;
    if (selector == "pilot") return std::string(test_case.id) == "development-strong-dense-0";
    if (selector == "development") return std::string(test_case.partition) == "development";
    if (selector == "heldout") return std::string(test_case.partition) == "heldout";
    const std::string prefix = "case:";
    return selector.rfind(prefix, 0) == 0 && std::string(test_case.id) == selector.substr(prefix.size());
}

std::vector<std::size_t> campaign_depths(const std::string& selector) {
    return selector == "pilot" ? std::vector<std::size_t>{32}
                                : std::vector<std::size_t>(DEPTHS.begin(), DEPTHS.end());
}

std::string active_name(const TestCase& test_case);
const char* regime_name(const TestCase& test_case);
std::vector<int> active_ties(const TestCase& test_case);
const char* topology_provenance(const TestCase& test_case);

void print_input_manifest(const TestCase& test_case) {
    double total = 0.0;
    for (double amount : test_case.amounts) total += amount;
    const auto fractions = composition(test_case.amounts);
    const auto sigma = epcsaft::segment_sigma(test_case.temperature, test_case.parameters);
    const auto diameter = epcsaft::segment_diameter(test_case.temperature, test_case.parameters, sigma);
    const double density = total / test_case.volume;
    const double eta = epcsaft::packing_fraction_from_geometry(density, fractions, test_case.parameters, diameter);
    const auto fields = recording_support::fields({
        recording_support::field("case_id", test_case.id),
        recording_support::field("partition", test_case.partition),
        recording_support::field("stratum", test_case.stratum),
        recording_support::field("regime", regime_name(test_case)),
        recording_support::field("topology_provenance", topology_provenance(test_case)),
        recording_support::field("k_ij", test_case.parameters.k_ij),
        recording_support::field("active_family", active_name(test_case)),
        recording_support::field("active_tied_indices", active_ties(test_case)),
        recording_support::field("generation_seed", test_case.seed.has_value() ? std::to_string(*test_case.seed) : "none"),
        recording_support::field("T_K", test_case.temperature),
        recording_support::field("V_m3", test_case.volume),
        recording_support::field("amounts_mol", test_case.amounts),
        recording_support::field("segment_count", test_case.parameters.segment_count),
        recording_support::field("sigma_constant_A", test_case.parameters.sigma_constant),
        recording_support::field("dispersion_energy_over_k_K", test_case.parameters.dispersion_energy_over_k),
        recording_support::field("site_component", test_case.parameters.site_component),
        recording_support::field("site_multiplicity", test_case.parameters.site_multiplicity),
        recording_support::field("association_active", test_case.parameters.association_active),
        recording_support::field("association_energy_over_k_K", test_case.parameters.association_energy_over_k),
        recording_support::field("association_kappa", test_case.parameters.association_volume),
        recording_support::field("sigma_A", sigma),
        recording_support::field("diameter_A", diameter),
        recording_support::field("packing_fraction", eta),
    });
    std::printf("EPCSAFT_ASSOCIATION_INPUT %s\n", recording_support::json_fields_object(fields).c_str());
    std::fflush(stdout);
}

void state_status(const TestCase& test_case, std::size_t depth, const char* status,
    const char* reason = nullptr) {
    std::printf("EPCSAFT_ASSOCIATION_STATE case_id=%s partition=%s stratum=%s depth=%zu status=%s",
        test_case.id, test_case.partition, test_case.stratum, depth, status);
    if (reason != nullptr) std::printf(" reason=%s", reason);
    std::putchar('\n');
    std::fflush(stdout);
}

std::vector<TestCase> test_cases() {
    auto weak = pure_two_b();
    weak.association_energy_over_k[1] = weak.association_energy_over_k[2] = 350.0;
    weak.association_volume[1] = weak.association_volume[2] = 0.01;
    auto strong = pure_two_b();
    strong.association_energy_over_k[1] = strong.association_energy_over_k[2] = 3200.0;
    strong.association_volume[1] = strong.association_volume[2] = 0.20;
    const auto associating = properties_support::associating_binary();
    std::vector<TestCase> result = {
        tuned("pure-2b", pure_two_b(), 300.0, 1.0e-4, {1.0}, std::nullopt, true),
        tuned("pure-3b", pure_three_b(), 300.0, 1.0e-4, {1.0}, std::nullopt, true),
        tuned("weak", weak, 450.0, 1.0e-4, {1.0}),
        tuned("strong", strong, 280.0, 1.0e-4, {1.0}),
        tuned("nonassociating-dilution", diluted_two_b(), 300.0, 1.0e-4, {0.05, 0.95}),
        tuned("trace-support", diluted_two_b(), 300.0, 1.0e-4, {1.0e-12, 1.0 - 1.0e-12}),
        tuned("disconnected", disconnected(), 315.0, 8.0e-5, {0.35, 0.65}),
        tuned("heterogeneous-cross", associating, 376.88, 1.0e-4, {0.1385, 0.8615}),
        tuned("active-energy", associating, 376.88, 1.0e-4, {0.1385, 0.8615},
            epcsaft::ActiveParameter{epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 4}}),
        tuned("active-volume", associating, 376.88, 1.0e-4, {0.1385, 0.8615},
            epcsaft::ActiveParameter{epcsaft::ParameterFamily::AssociationVolume, {1, 4}}),
        tuned("neutral-three-component", neutral_three_component(), 325.0, 7.5e-5, {0.2, 0.5, 0.3}),
        seeded_case("seeded-mixture-a", 21701),
        seeded_case("seeded-mixture-b", 21702),
        tuned("exact-zero-support", diluted_two_b(), 300.0, 1.0e-4, {0.0, 1.0}),
    };
    const auto generated = generated_cases();
    result.insert(result.end(), generated.begin(), generated.end());
    return result;
}

std::string active_name(const TestCase& test_case) {
    if (!test_case.active.has_value()) return "none";
    return test_case.active->family == epcsaft::ParameterFamily::AssociationEnergyOverK
        ? "association_energy_over_k" : "association_volume";
}

const char* regime_name(const TestCase& test_case) {
    const std::string stratum = test_case.stratum;
    if (stratum == "trace-support" || stratum == "exact-zero-support") return "support";
    if (stratum == "parameter-boundaries") return "boundary";
    if (stratum == "strong-dense" || stratum == "heterogeneous-cross-capacity"
        || stratum == "disconnected-nonassociating") return "stress";
    return stratum == "controls" ? "control" : "ordinary";
}

std::vector<int> active_ties(const TestCase& test_case) { if (!test_case.active.has_value()) return {}; std::vector<int> result; for (const std::size_t index : test_case.active->indices) result.push_back(static_cast<int>(index)); return result; }

const char* topology_provenance(const TestCase& test_case) {
    const std::string id = test_case.id;
    if (test_case.seed.has_value()) return "manufactured-seeded";
    return id == "heterogeneous-cross" || id == "active-energy" || id == "active-volume"
        ? "source-fixture" : "manufactured";
}

std::string inputs_for(const TestCase& test_case, const Candidate&, const Reference& reference,
    std::size_t depth) {
    return recording_support::fields({recording_support::field("case_id", test_case.id),
        recording_support::field("partition", test_case.partition), recording_support::field("stratum", test_case.stratum),
        recording_support::field("regime", regime_name(test_case)),
        recording_support::field("depth", depth),
        recording_support::field("T_K", test_case.temperature), recording_support::field("V_m3", test_case.volume),
        recording_support::field("amounts_mol", test_case.amounts), recording_support::field("active_parameter", active_name(test_case)),
        recording_support::field("generation_seed", test_case.seed.has_value() ? std::to_string(*test_case.seed) : "none"),
        recording_support::field("reference_value_status", reference.value_status),
        recording_support::field("reference_derivative_status", reference.derivative_status),
        recording_support::field("reference_error", reference.error),
        recording_support::field("reference_newton_defect", reference.newton_defect)});
}

double scale_for(const TestCase& test_case, std::size_t coordinate) {
    const std::size_t components = test_case.parameters.component_count();
    if (coordinate == 0) return test_case.temperature;
    if (coordinate == 1) return test_case.volume;
    if (coordinate < 2 + components) {
        double total = 0.0;
        for (double amount : test_case.amounts) total += amount;
        return total;
    }
    return test_case.active->family == epcsaft::ParameterFamily::AssociationEnergyOverK ? 1000.0 : 0.1;
}

void reference_diagnostics(const TestCase& test_case, const Reference& reference, const std::string& inputs) {
    const auto diagnostic_inputs = recording_support::fields({inputs,
        recording_support::field("worst_first_coordinate", reference.association_tensor.first_sensitivity_coordinate),
        recording_support::field("worst_second_row", reference.association_tensor.second_sensitivity_row),
        recording_support::field("worst_second_column", reference.association_tensor.second_sensitivity_column)});
    recording_support::observe((std::string("reference-quality/") + test_case.id).c_str(), "reference_condition", "1",
        reference.condition, reference.condition, diagnostic_inputs, "Diagnostic", "implicit_reference_quality", reference.status.c_str(),
        "reference-quality", test_case.id, "condition", "T,V,n,active", "cond_inf(J*diag(X))", "case");
    const auto& tensor = reference.association_tensor;
    recording_support::observe((std::string("reference-first-sensitivity/") + test_case.id).c_str(),
        "first_sensitivity_residual", "1", tensor.first_sensitivity_residual, tensor.first_sensitivity_residual,
        diagnostic_inputs, "Diagnostic", "implicit_reference_quality", reference.derivative_status.c_str(),
        "reference-first-sensitivity", test_case.id, "R_1", "T,V,n,active", "scaled JX_q+F_q", "coordinate");
    recording_support::observe((std::string("reference-second-sensitivity/") + test_case.id).c_str(),
        "second_sensitivity_residual", "1", tensor.second_sensitivity_residual, tensor.second_sensitivity_residual,
        diagnostic_inputs, "Diagnostic", "implicit_reference_quality", reference.derivative_status.c_str(),
        "reference-second-sensitivity", test_case.id, "R_2", "T,V,n,active", "scaled JX_qr+contraction", "row,column");
    recording_support::observe((std::string("reference-sensitivity-proxy/") + test_case.id).c_str(),
        "sensitivity_error_proxy", "1", tensor.sensitivity_error_proxy, tensor.sensitivity_error_proxy,
        diagnostic_inputs, "Diagnostic", "implicit_reference_quality", reference.derivative_status.c_str(),
        "reference-sensitivity-proxy", test_case.id, "kappa_inf*(R+n eps)", "T,V,n,active", "screening proxy", "case");
    recording_support::observe((std::string("reference-newton-defect/") + test_case.id).c_str(),
        "newton_defect", "1", reference.newton_defect, reference.newton_defect,
        diagnostic_inputs, "Diagnostic", "implicit_reference_quality", reference.value_status.c_str(),
        "reference-newton-defect", test_case.id, "D_X^-1 J^-1 F", "T,V,n,active", "site-value screen", "case");
}

bool compare_row(const char* comparison_id, const char* quantity, const char* unit, double model, double reference,
    const std::string& inputs, const char* basis, const char* method, const char* claim, const char* criterion,
    double atol, double rtol, const char* comparison_group_id = nullptr, const char* point_id = nullptr,
    const char* quantity_latex = nullptr, const char* inputs_latex = nullptr, const char* relationship_latex = nullptr,
    const char* coordinates = nullptr) {
    const bool verdict = recording_support::scalar_verdict(model, reference, criterion, atol, rtol);
    recording_support::print_json_compare(comparison_id, quantity, unit, model, reference, inputs, basis, method, claim,
        criterion, atol, rtol, verdict, verdict ? "Criterion satisfied" : "Criterion failed", comparison_group_id,
        point_id, quantity_latex, inputs_latex, relationship_latex, coordinates);
    return verdict;
}

std::string depth_id(std::size_t depth, const char* group, const std::string& point) {
    return "d" + std::to_string(depth) + "/" + group + "/" + point;
}

std::string depth_prefix(std::size_t depth, const char* group) {
    return "d" + std::to_string(depth) + "/" + group + "/";
}

bool compare_values(const TestCase& test_case, const Candidate& candidate, const Reference& reference,
    std::size_t depth) {
    bool passed = true;
    const auto inputs = inputs_for(test_case, candidate, reference, depth);
    for (std::size_t site = 0; site < candidate.sites.size(); ++site) {
        const std::string point = std::string(test_case.id) + "/" + std::to_string(site);
        const std::string id = depth_id(depth, "site-fraction", point);
        passed = compare_row(id.c_str(), "site_fraction", "1",
            candidate.sites[site], reference.sites[site], inputs, "Comparison", "fixed_depth_candidate", "site fraction agreement",
            "reference_scale", 1.0e-12, SITE_RTOL, "site-fraction", ("d" + std::to_string(depth) + "/" + point).c_str(), "X", "T,V,n", "X_candidate=X_reference", "site") && passed;
        const std::string residual_id = depth_id(depth, "candidate-residual", point);
        passed = compare_row(residual_id.c_str(), "mass_action_residual", "1",
            candidate.mass_action[site], 0.0, inputs, "Comparison", "defining_residual", "candidate original residual",
            "strict_absolute", VALUE_TOL, 0.0, "candidate-residual", ("d" + std::to_string(depth) + "/" + point).c_str(), "F", "T,V,n", "F_candidate=0", "site") && passed;
        const std::string reference_id = depth_id(depth, "reference-residual", point);
        passed = compare_row(reference_id.c_str(), "reference_mass_action_residual", "1",
            reference.mass_action[site], 0.0, inputs, "Reference", "defining_residual", "reference original residual",
            "strict_absolute", REFERENCE_TOL, 0.0, "reference-residual", ("d" + std::to_string(depth) + "/" + point).c_str(), "F", "T,V,n", "F_reference=0", "site") && passed;
    }
    const std::string association_id = depth_id(depth, "association-energy", test_case.id);
    passed = compare_row(association_id.c_str(), "Phi_association", "1",
        candidate.association, reference.association, inputs, "Comparison", "fixed_depth_candidate", "extensive association energy",
        "strict_reference_or_one", 0.0, ENERGY_RTOL, "association-energy", ("d" + std::to_string(depth) + "/" + test_case.id).c_str(), "Phi_assoc", "T,V,n", "candidate=reference", "case") && passed;
    const std::string residual_id = depth_id(depth, "residual-energy", test_case.id);
    passed = compare_row(residual_id.c_str(), "Phi_residual", "1",
        candidate.residual, reference.residual, inputs, "Comparison", "fixed_depth_candidate", "extensive residual energy",
        "strict_reference_or_one", 0.0, ENERGY_RTOL, "residual-energy", ("d" + std::to_string(depth) + "/" + test_case.id).c_str(), "Phi_res", "T,V,n", "candidate=reference", "case") && passed;
    return passed;
}

bool compare_derivative_set(const TestCase& test_case, const Candidate& candidate, const Reference& reference,
    const std::vector<double>& model_gradient, const std::vector<double>& model_hessian,
    const std::vector<double>& reference_gradient, const std::vector<double>& reference_hessian,
    const char* label, const char* claim, std::size_t depth) {
    const std::size_t width = model_gradient.size();
    const auto inputs = inputs_for(test_case, candidate, reference, depth);
    bool passed = true;
    const std::string gradient_group = std::string(label) + "-gradient";
    const std::string hessian_group = std::string(label) + "-hessian";
    for (std::size_t coordinate = 0; coordinate < width; ++coordinate) {
        const double scale = scale_for(test_case, coordinate);
        const std::string point = std::string(test_case.id) + "/" + std::to_string(coordinate);
        const auto details = recording_support::fields({inputs, recording_support::field("coordinate", coordinate),
            recording_support::field("scale", scale), recording_support::field("model_raw", model_gradient[coordinate]),
            recording_support::field("reference_raw", reference_gradient[coordinate])});
        const std::string id = depth_id(depth, gradient_group.c_str(), point);
        passed = compare_row(id.c_str(), "scaled_gradient", "1",
            model_gradient[coordinate] * scale, reference_gradient[coordinate] * scale, details, "Comparison", "CppAD_vs_implicit_reference",
            claim, "reference_scale", DERIVATIVE_ATOL, DERIVATIVE_RTOL, gradient_group.c_str(), ("d" + std::to_string(depth) + "/" + point).c_str(), "dPhi/dq_scaled", "T,V,n,active", "candidate=reference", "coordinate") && passed;
    }
    for (std::size_t row = 0; row < width; ++row) for (std::size_t column = 0; column < width; ++column) {
        const std::size_t entry = row * width + column;
        const double scale = scale_for(test_case, row) * scale_for(test_case, column);
        const std::string point = std::string(test_case.id) + "/" + std::to_string(row) + "-" + std::to_string(column);
        const auto details = recording_support::fields({inputs, recording_support::field("row", row), recording_support::field("column", column),
            recording_support::field("scale", scale), recording_support::field("model_raw", model_hessian[entry]),
            recording_support::field("reference_raw", reference_hessian[entry])});
            const std::string id = depth_id(depth, hessian_group.c_str(), point);
            passed = compare_row(id.c_str(), "scaled_hessian", "1",
            model_hessian[entry] * scale, reference_hessian[entry] * scale, details, "Comparison", "CppAD_vs_implicit_reference",
            claim, "reference_scale", DERIVATIVE_ATOL, DERIVATIVE_RTOL, hessian_group.c_str(), ("d" + std::to_string(depth) + "/" + point).c_str(), "d2Phi/dq_scaled2", "T,V,n,active", "candidate=reference", "row,column") && passed;
    }
    return passed;
}

EnergyValues candidate_energy(const TestCase& test_case, std::size_t depth = QUALIFYING_DEPTH) {
    try {
        Tape tape = make_tape(test_case, depth);
        const auto values = tape.function.Forward(0, tape.point);
        return {values[0], values[1], std::isfinite(values[0]) && std::isfinite(values[1])};
    } catch (...) {
        AD::abort_recording();
        return {};
    }
}

EnergyValues reference_energy(const TestCase& test_case) {
    try {
        const Reference reference = evaluate_reference(test_case, false);
        return {reference.association, reference.residual, reference.value_available
            && std::isfinite(reference.association) && std::isfinite(reference.residual)};
    } catch (...) {
        AD::abort_recording();
        return {};
    }
}

TestCase shifted_case(const TestCase& test_case, std::size_t coordinate, double amount) {
    TestCase result = test_case;
    if (coordinate == 0) result.temperature += amount;
    else if (coordinate == 1) result.volume += amount;
    else if (coordinate < 2 + result.parameters.component_count()) {
        result.amounts[coordinate - 2] += amount;
    } else {
        auto& values = result.active->family == epcsaft::ParameterFamily::AssociationEnergyOverK
            ? result.parameters.association_energy_over_k : result.parameters.association_volume;
        for (const std::size_t index : result.active->indices) values[index] += amount;
    }
    return result;
}

bool valid_shift(const TestCase& test_case, std::size_t coordinate, double amount) {
    const TestCase shifted = shifted_case(test_case, coordinate, amount);
    if (!(shifted.temperature > 0.0) || !(shifted.volume > 0.0)) return false;
    double total = 0.0;
    for (const double value : shifted.amounts) {
        if (!std::isfinite(value) || value < 0.0) return false;
        total += value;
    }
    if (!(total > 0.0)) return false;
    if (shifted.active.has_value()) {
        const auto& values = shifted.active->family == epcsaft::ParameterFamily::AssociationEnergyOverK
            ? shifted.parameters.association_energy_over_k : shifted.parameters.association_volume;
        for (const std::size_t index : shifted.active->indices)
            if (!std::isfinite(values[index]) || values[index] < 0.0) return false;
    }
    return true;
}

struct Difference final {
    double first = std::numeric_limits<double>::quiet_NaN();
    double second = std::numeric_limits<double>::quiet_NaN();
    const char* stencil = "unavailable";
    bool available = false;
};

using EnergyEvaluator = EnergyValues (*)(const TestCase&, std::size_t);

Difference finite_difference(EnergyEvaluator evaluator, const TestCase& test_case,
    std::size_t coordinate, double step, int energy) {
    const EnergyValues center_values = evaluator(test_case, QUALIFYING_DEPTH);
    if (!center_values.available) return {};
    bool available = true;
    const auto value = [&](double offset) {
        const EnergyValues values = evaluator(shifted_case(test_case, coordinate, offset), QUALIFYING_DEPTH);
        available = available && values.available;
        return energy == 0 ? values.association : values.residual;
    };
    const bool central = valid_shift(test_case, coordinate, -2.0 * step)
        && valid_shift(test_case, coordinate, 2.0 * step);
    const bool forward = valid_shift(test_case, coordinate, 4.0 * step);
    const bool backward = valid_shift(test_case, coordinate, -4.0 * step);
    Difference result;
    const double center = energy == 0 ? center_values.association : center_values.residual;
    if (central) {
        const double m2 = value(-2.0 * step), m1 = value(-step);
        const double p1 = value(step), p2 = value(2.0 * step);
        result.first = (m2 - 8.0 * m1 + 8.0 * p1 - p2) / (12.0 * step);
        result.second = (-m2 + 16.0 * m1 - 30.0 * center + 16.0 * p1 - p2) / (12.0 * step * step);
        result.stencil = "five_point_central";
    } else if (forward) {
        const double p1 = value(step), p2 = value(2.0 * step);
        const double p3 = value(3.0 * step), p4 = value(4.0 * step);
        result.first = (-25.0 * center + 48.0 * p1 - 36.0 * p2 + 16.0 * p3 - 3.0 * p4) / (12.0 * step);
        result.second = (35.0 * center - 104.0 * p1 + 114.0 * p2 - 56.0 * p3 + 11.0 * p4) / (12.0 * step * step);
        result.stencil = "five_point_forward";
    } else if (backward) {
        const double m1 = value(-step), m2 = value(-2.0 * step);
        const double m3 = value(-3.0 * step), m4 = value(-4.0 * step);
        result.first = (25.0 * center - 48.0 * m1 + 36.0 * m2 - 16.0 * m3 + 3.0 * m4) / (12.0 * step);
        result.second = (35.0 * center - 104.0 * m1 + 114.0 * m2 - 56.0 * m3 + 11.0 * m4) / (12.0 * step * step);
        result.stencil = "five_point_backward";
    } else return result;
    result.available = available && std::isfinite(result.first) && std::isfinite(result.second);
    return result;
}

Difference mixed_difference(EnergyEvaluator evaluator, const TestCase& test_case,
    std::size_t row, std::size_t column, double row_step, double column_step, int energy) {
    const EnergyValues center_values = evaluator(test_case, QUALIFYING_DEPTH);
    if (!center_values.available) return {};
    bool available = true;
    const auto value = [&](double row_offset, double column_offset) {
        TestCase shifted = shifted_case(test_case, row, row_offset);
        if (column == 0) shifted.temperature += column_offset;
        else if (column == 1) shifted.volume += column_offset;
        else if (column < 2 + shifted.parameters.component_count()) shifted.amounts[column - 2] += column_offset;
        else {
            auto& values = shifted.active->family == epcsaft::ParameterFamily::AssociationEnergyOverK
                ? shifted.parameters.association_energy_over_k : shifted.parameters.association_volume;
            for (const std::size_t index : shifted.active->indices) values[index] += column_offset;
        }
        const EnergyValues values = evaluator(shifted, QUALIFYING_DEPTH);
        available = available && values.available;
        return energy == 0 ? values.association : values.residual;
    };
    const auto valid_pair = [&](double row_offset, double column_offset) {
        return valid_shift(test_case, row, row_offset)
            && valid_shift(shifted_case(test_case, row, row_offset), column, column_offset);
    };
    Difference result;
    if (!valid_pair(row_step, column_step) || !valid_pair(row_step, -column_step)
        || !valid_pair(-row_step, column_step) || !valid_pair(-row_step, -column_step)) return result;
    result.first = (value(row_step, column_step) - value(row_step, -column_step)
        - value(-row_step, column_step) + value(-row_step, -column_step))
        / (4.0 * row_step * column_step);
    result.available = available && std::isfinite(result.first);
    result.stencil = "four_point_central";
    return result;
}

void observe_centered_unavailable(const std::string& id, const char* quantity,
    const std::string& inputs, const char* group, const char* point,
    const char* relationship, const char* coordinates) {
    recording_support::observe(id.c_str(), quantity, "1", std::numeric_limits<double>::quiet_NaN(),
        std::nullopt, inputs, "Independent_check", "finite_difference", "finite-difference check unavailable",
        group, point, quantity[9] == 'g' ? "dPhi/dq_scaled" : "d2Phi/dq_scaled2", "T,V,n,active", relationship, coordinates);
}

bool centered_checks(const TestCase& test_case, const std::vector<std::size_t>& coordinates,
    const std::vector<std::pair<std::size_t, std::size_t>>& mixed_pairs,
    const std::array<double, 2>& relative_steps = {1.0e-4, 3.0e-5}) {
    const std::size_t width = 2 + test_case.parameters.component_count()
        + (test_case.active.has_value() ? 1 : 0);
    Candidate candidate;
    Reference reference;
    bool center_candidate = true, center_reference = true;
    try { candidate = evaluate_candidate(test_case, QUALIFYING_DEPTH); }
    catch (...) { AD::abort_recording(); center_candidate = false; }
    try { reference = evaluate_reference(test_case); }
    catch (...) { AD::abort_recording(); center_reference = false; }
    const bool candidate_derivatives = center_candidate && candidate.association_gradient.size() == width
        && candidate.residual_gradient.size() == width && candidate.association_hessian.size() == width * width
        && candidate.residual_hessian.size() == width * width;
    const bool reference_derivatives = center_reference && reference.derivative_available
        && reference.association_tensor.gradient.size() == width && reference.tensor.gradient.size() == width
        && reference.association_tensor.hessian.size() == width * width && reference.tensor.hessian.size() == width * width;
    bool passed = candidate_derivatives && reference_derivatives;
    for (const double relative_step : relative_steps) {
        for (const std::size_t coordinate : coordinates) {
            const double scale = scale_for(test_case, coordinate);
            const double step = relative_step * scale;
            for (int energy = 0; energy < 2; ++energy) {
                const char* energy_name = energy == 0 ? "association" : "residual";
                const auto inputs = recording_support::fields({recording_support::field("case_id", test_case.id),
                    recording_support::field("source", "candidate_and_reference"), recording_support::field("energy", energy_name),
                    recording_support::field("coordinate", coordinate), recording_support::field("scale", scale),
                    recording_support::field("step", step), recording_support::field("relative_step", relative_step)});
                const std::string point = std::string(test_case.id) + "/" + energy_name + "/h"
                    + recording_support::number(relative_step) + "/" + std::to_string(coordinate);
                const Difference candidate_difference = finite_difference(
                    candidate_energy, test_case, coordinate, step, energy);
                const Difference reference_difference = finite_difference(
                    [](const TestCase& shifted, std::size_t) { return reference_energy(shifted); },
                    test_case, coordinate, step, energy);
                const std::vector<double>& candidate_gradient = energy == 0 ? candidate.association_gradient : candidate.residual_gradient;
                const std::vector<double>& candidate_hessian = energy == 0 ? candidate.association_hessian : candidate.residual_hessian;
                const std::vector<double>& reference_gradient = energy == 0 ? reference.association_tensor.gradient : reference.tensor.gradient;
                const std::vector<double>& reference_hessian = energy == 0 ? reference.association_tensor.hessian : reference.tensor.hessian;
                const auto compare_difference = [&](const Difference& difference, double analytic_first, double analytic_second,
                    const char* source, bool source_available) {
                    const auto row_inputs = recording_support::fields({inputs, recording_support::field("source", source),
                        recording_support::field("stencil", difference.stencil), recording_support::field("status", difference.available && source_available ? "compared" : "unavailable")});
                    const std::string prefix = std::string("centered-gradient/") + source + "/" + point;
                    if (difference.available && source_available) {
                        bool result = compare_row(prefix.c_str(), "centered_gradient", "1", difference.first * scale,
                            analytic_first * scale, row_inputs, "Independent_check", "finite_difference", "finite-difference energy gradient",
                            "reference_scale", CENTERED_ATOL, CENTERED_RTOL, "centered-gradient", point.c_str(), "dPhi/dq_scaled", "T,V,n,active", "stencil=graph", "coordinate");
                        const std::string hessian_id = std::string("centered-hessian/") + source + "/" + point;
                        result = compare_row(hessian_id.c_str(), "centered_hessian", "1", difference.second * scale * scale,
                            analytic_second * scale * scale, row_inputs, "Independent_check", "finite_difference", "finite-difference energy Hessian",
                            "reference_scale", CENTERED_ATOL, CENTERED_RTOL, "centered-hessian", point.c_str(), "d2Phi/dq_scaled2", "T,V,n,active", "stencil=graph", "coordinate") && result;
                        return result;
                    }
                    observe_centered_unavailable(prefix, "centered_gradient", row_inputs, "centered-gradient", point.c_str(), "finite_difference", "coordinate");
                    observe_centered_unavailable(std::string("centered-hessian/") + source + "/" + point, "centered_hessian", row_inputs, "centered-hessian", point.c_str(), "finite_difference", "coordinate");
                    return false;
                };
                passed = compare_difference(candidate_difference,
                    candidate_derivatives ? candidate_gradient[coordinate] : std::numeric_limits<double>::quiet_NaN(),
                    candidate_derivatives ? candidate_hessian[coordinate * width + coordinate] : std::numeric_limits<double>::quiet_NaN(),
                    "candidate", candidate_derivatives) && passed;
                passed = compare_difference(reference_difference,
                    reference_derivatives ? reference_gradient[coordinate] : std::numeric_limits<double>::quiet_NaN(),
                    reference_derivatives ? reference_hessian[coordinate * width + coordinate] : std::numeric_limits<double>::quiet_NaN(),
                    "reference", reference_derivatives) && passed;
            }
        }
        for (const auto [row, column] : mixed_pairs) for (int energy = 0; energy < 2; ++energy) {
            const char* energy_name = energy == 0 ? "association" : "residual";
            const double row_step = relative_step * scale_for(test_case, row);
            const double column_step = relative_step * scale_for(test_case, column);
            const auto inputs = recording_support::fields({recording_support::field("case_id", test_case.id),
                recording_support::field("energy", energy_name), recording_support::field("row", row), recording_support::field("column", column),
                recording_support::field("row_step", row_step), recording_support::field("column_step", column_step),
                recording_support::field("relative_step", relative_step)});
            const std::string point = std::string(test_case.id) + "/" + energy_name + "/h"
                + recording_support::number(relative_step) + "/" + std::to_string(row) + "-" + std::to_string(column);
            const Difference candidate_difference = mixed_difference(candidate_energy, test_case, row, column, row_step, column_step, energy);
            const Difference reference_difference = mixed_difference(
                [](const TestCase& shifted, std::size_t) { return reference_energy(shifted); },
                test_case, row, column, row_step, column_step, energy);
            const auto compare_mixed = [&](const Difference& difference, double analytic, const char* source, bool source_available) {
                const auto row_inputs = recording_support::fields({inputs, recording_support::field("source", source),
                    recording_support::field("stencil", difference.stencil), recording_support::field("status", difference.available && source_available ? "compared" : "unavailable")});
                const std::string id = std::string("centered-mixed/") + source + "/" + point;
                if (difference.available && source_available)
                    return compare_row(id.c_str(), "centered_mixed_hessian", "1",
                        difference.first * scale_for(test_case, row) * scale_for(test_case, column),
                        analytic * scale_for(test_case, row) * scale_for(test_case, column), row_inputs,
                        "Independent_check", "finite_difference", "finite-difference mixed energy Hessian", "reference_scale",
                        CENTERED_ATOL, CENTERED_RTOL, "centered-mixed", point.c_str(), "d2Phi/dq_scaled2", "T,V,n,active", "stencil=graph", "row,column");
                observe_centered_unavailable(id, "centered_mixed_hessian", row_inputs, "centered-mixed", point.c_str(), "finite_difference", "row,column");
                return false;
            };
            const std::vector<double>& candidate_hessian = energy == 0 ? candidate.association_hessian : candidate.residual_hessian;
            const std::vector<double>& reference_hessian = energy == 0 ? reference.association_tensor.hessian : reference.tensor.hessian;
            passed = compare_mixed(candidate_difference, candidate_derivatives ? candidate_hessian[row * width + column] : std::numeric_limits<double>::quiet_NaN(), "candidate", candidate_derivatives) && passed;
            passed = compare_mixed(reference_difference, reference_derivatives ? reference_hessian[row * width + column] : std::numeric_limits<double>::quiet_NaN(), "reference", reference_derivatives) && passed;
        }
    }
    return passed;
}

bool analytic_controls(const TestCase& test_case, const Candidate& candidate, bool candidate_ok,
    std::size_t depth) {
    if (!test_case.analytic) return true;
    const std::size_t expected_size = test_case.parameters.site_count();
    if (!candidate_ok || candidate.sites.size() != expected_size) {
        for (std::size_t site = 0; site < expected_size; ++site)
            recording_support::observe(depth_id(depth, "analytic-control", std::string(test_case.id) + "/" + std::to_string(site)).c_str(),
                "site_fraction", "1", std::numeric_limits<double>::quiet_NaN(), std::nullopt,
                recording_support::fields({recording_support::field("case_id", test_case.id), recording_support::field("status", "candidate_unavailable")}),
                "Analytic_control", "closed_form", "candidate matches pure-site control", "analytic-control", test_case.id, "X", "T,V,n", "status", "site");
        return false;
    }
    const auto mole_fractions = composition(test_case.amounts);
    const double density = test_case.amounts.front() / test_case.volume;
    const auto sigma = epcsaft::segment_sigma(test_case.temperature, test_case.parameters);
    const auto diameter = epcsaft::segment_diameter(test_case.temperature, test_case.parameters, sigma);
    const auto delta = epcsaft::association_delta(test_case.temperature, density, mole_fractions, test_case.parameters, sigma, diameter);
    const auto k = association_k(density, mole_fractions, test_case.parameters, delta);
    std::vector<double> expected(candidate.sites.size());
    if (expected.size() == 2) {
        const double root = 2.0 / (1.0 + std::sqrt(1.0 + 4.0 * k[1]));
        expected = {root, root};
    } else {
        const double donor = 2.0 / (1.0 + k[1] + std::sqrt(1.0 + 6.0 * k[1] + k[1] * k[1]));
        expected = {donor, (1.0 + donor) / 2.0, (1.0 + donor) / 2.0};
    }
    bool passed = true;
    for (std::size_t site = 0; site < expected.size(); ++site) {
        const std::string id = depth_id(depth, "analytic-control", std::string(test_case.id) + "/" + std::to_string(site));
        passed = compare_row(id.c_str(),
            "site_fraction", "1", candidate.sites[site], expected[site], recording_support::fields({recording_support::field("case_id", test_case.id), recording_support::field("coupling", k[1])}),
            "Analytic_control", "closed_form", "candidate matches pure-site control", "strict_reference_scale", 1.0e-12, 1.0e-12,
            "analytic-control", test_case.id, "X", "T,V,n", "candidate=closed_form", "site") && passed;
    }
    return passed;
}

bool k_identity(const TestCase& test_case, const Candidate& candidate, std::size_t depth) {
    const auto mole_fractions = composition(test_case.amounts);
    double total = 0.0;
    for (double amount : test_case.amounts) total += amount;
    const double density = total / test_case.volume;
    const auto sigma = epcsaft::segment_sigma(test_case.temperature, test_case.parameters);
    const auto diameter = epcsaft::segment_diameter(test_case.temperature, test_case.parameters, sigma);
    const auto delta = epcsaft::association_delta(test_case.temperature, density, mole_fractions, test_case.parameters, sigma, diameter);
    const auto k = association_k(density, mole_fractions, test_case.parameters, delta);
    bool passed = true;
    for (std::size_t row = 0; row < candidate.sites.size(); ++row) {
        double bonded = 0.0;
        for (std::size_t column = 0; column < candidate.sites.size(); ++column)
            bonded += k[row * candidate.sites.size() + column] * candidate.sites[column];
        const double manual = candidate.sites[row] * (1.0 + bonded) - 1.0;
        const std::string id = depth_id(depth, "k-residual-identity", std::string(test_case.id) + "/" + std::to_string(row));
        passed = compare_row(id.c_str(),
            "mass_action_residual", "1", candidate.mass_action[row], manual,
            recording_support::fields({recording_support::field("case_id", test_case.id), recording_support::field("K_row", std::vector<double>(k.begin() + row * candidate.sites.size(), k.begin() + (row + 1) * candidate.sites.size()))}),
            "Relationship", "shared_residual", "association residual uses rho_hat x_b nu_b Delta", "strict_reference_scale", 1.0e-14, 1.0e-14,
            "k-residual-identity", test_case.id, "F", "T,V,n", "F_shared=F_manual", "site") && passed;
    }
    return passed;
}

void emit_derivative_unavailable(const TestCase& test_case, const Candidate& candidate,
    const Reference& reference, std::size_t depth, const char* status);

void emit_unavailable(const TestCase& test_case, const Candidate& candidate, const Reference& reference,
    std::size_t depth) {
    const auto inputs = inputs_for(test_case, candidate, reference, depth);
    const std::size_t sites = test_case.parameters.site_count();
    for (std::size_t site = 0; site < sites; ++site) {
        const auto candidate_site = site < candidate.sites.size() ? std::optional<double>(candidate.sites[site]) : std::nullopt;
        const auto reference_site = site < reference.sites.size() ? std::optional<double>(reference.sites[site]) : std::nullopt;
        recording_support::observe(depth_id(depth, "site-fraction", std::string(test_case.id) + "/" + std::to_string(site)).c_str(), "site_fraction", "1",
            candidate_site.value_or(std::numeric_limits<double>::quiet_NaN()), reference_site, inputs, "Comparison", "fixed_depth_candidate", reference.status.c_str(), "site-fraction", test_case.id, "X", "T,V,n", "status", "site");
        recording_support::observe(depth_id(depth, "candidate-residual", std::string(test_case.id) + "/" + std::to_string(site)).c_str(), "mass_action_residual", "1",
            site < candidate.mass_action.size() ? candidate.mass_action[site] : std::numeric_limits<double>::quiet_NaN(), std::nullopt, inputs, "Comparison", "defining_residual", reference.status.c_str(), "candidate-residual", test_case.id, "F", "T,V,n", "status", "site");
        recording_support::observe(depth_id(depth, "reference-residual", std::string(test_case.id) + "/" + std::to_string(site)).c_str(), "reference_mass_action_residual", "1",
            site < reference.mass_action.size() ? reference.mass_action[site] : std::numeric_limits<double>::quiet_NaN(), std::nullopt, inputs, "Reference", "defining_residual", reference.status.c_str(), "reference-residual", test_case.id, "F", "T,V,n", "status", "site");
    }
    recording_support::observe(depth_id(depth, "association-energy", test_case.id).c_str(), "Phi_association", "1", candidate.association, std::nullopt, inputs, "Comparison", "fixed_depth_candidate", reference.status.c_str(), "association-energy", test_case.id, "Phi_assoc", "T,V,n", "status", "case");
    recording_support::observe(depth_id(depth, "residual-energy", test_case.id).c_str(), "Phi_residual", "1", candidate.residual, std::nullopt, inputs, "Comparison", "fixed_depth_candidate", reference.status.c_str(), "residual-energy", test_case.id, "Phi_res", "T,V,n", "status", "case");
    emit_derivative_unavailable(test_case, candidate, reference, depth, reference.status.c_str());
    for (std::size_t site = 0; site < sites; ++site)
        recording_support::observe(depth_id(depth, "k-residual-identity", std::string(test_case.id) + "/" + std::to_string(site)).c_str(),
            "mass_action_residual", "1", std::numeric_limits<double>::quiet_NaN(), std::nullopt, inputs,
            "Relationship", "shared_residual", reference.status.c_str(), "k-residual-identity", test_case.id, "F", "T,V,n", "status", "site");
}

void emit_derivative_unavailable(const TestCase& test_case, const Candidate& candidate,
    const Reference& reference, std::size_t depth, const char* status) {
    const auto inputs = inputs_for(test_case, candidate, reference, depth);
    const std::size_t width = 2 + test_case.parameters.component_count()
        + (test_case.active.has_value() ? 1 : 0);
    for (const char* label : {"association", "residual"}) {
        const std::string gradient_group = std::string(label) + "-gradient";
        const std::string hessian_group = std::string(label) + "-hessian";
        for (std::size_t coordinate = 0; coordinate < width; ++coordinate)
            recording_support::observe(depth_id(depth, gradient_group.c_str(),
                std::string(test_case.id) + "/" + std::to_string(coordinate)).c_str(),
                "scaled_gradient", "1", std::numeric_limits<double>::quiet_NaN(), std::nullopt,
                recording_support::fields({inputs, recording_support::field("status", status)}),
                "Comparison", "CppAD_vs_implicit_reference", status, gradient_group.c_str(), test_case.id,
                "dPhi/dq", "T,V,n,active", "status", "coordinate");
        for (std::size_t entry = 0; entry < width * width; ++entry)
            recording_support::observe(depth_id(depth, hessian_group.c_str(),
                std::string(test_case.id) + "/" + std::to_string(entry)).c_str(),
                "scaled_hessian", "1", std::numeric_limits<double>::quiet_NaN(), std::nullopt,
                recording_support::fields({inputs, recording_support::field("status", status)}),
                "Comparison", "CppAD_vs_implicit_reference", status, hessian_group.c_str(), test_case.id,
                "d2Phi/dq2", "T,V,n,active", "status", "entry");
    }
}

bool production_zero_support_status(const TestCase& test_case) {
    if (std::string(test_case.id) != "exact-zero-support") return true;
    // The former singular-Hessian limitation at an exactly absent associating species is
    // resolved in production; the probe now requires a finite Hessian there.
    bool finite = false;
    std::string error_message;
    try {
        epcsaft::HelmholtzModel model(test_case.parameters);
        const auto result = model.evaluate(test_case.temperature, test_case.volume, test_case.amounts, epcsaft::DerivativeOrder::Hessian);
        finite = std::all_of(result.hessian.begin(), result.hessian.end(), [](double v) { return std::isfinite(v); });
    } catch (const std::exception& error) {
        error_message = error.what();
    }
    recording_support::observe("production-weighted-hessian/exact-zero-support", "status", "1",
        finite ? 1.0 : 0.0, std::nullopt, "case_id:exact-zero-support;exception=" + error_message,
        "Diagnostic", "production_probe", "exact-zero support Hessian is finite", "production-weighted-hessian", "exact-zero-support", "status", "T,V,n", "finite Hessian", "case");
    return finite;
}

bool extensive_scaling(const TestCase& test_case) {
    TestCase scaled = test_case;
    for (double& amount : scaled.amounts) amount *= 2.0;
    scaled.volume *= 2.0;
    const double center = candidate_energy(test_case).residual;
    const double doubled = candidate_energy(scaled).residual;
    return compare_row("extensive-scaling/heterogeneous-cross", "Phi_residual", "1", doubled, 2.0 * center,
        "case_id:heterogeneous-cross;amount_scale:2;volume_scale:2", "Relationship", "candidate_extensivity", "Phi_res scales with total amount at fixed density",
        "strict_reference_scale", 1.0e-10, 1.0e-12, "extensive-scaling", "heterogeneous-cross", "Phi_res", "T,V,n", "Phi(2n,2V)=2Phi(n,V)", "case");
}

bool failure_continuation(const TestCase& valid_case) {
    TestCase invalid = valid_case;
    std::fill(invalid.amounts.begin(), invalid.amounts.end(), 0.0);
    bool failed = false;
    try {
        const Candidate candidate = evaluate_candidate(invalid, QUALIFYING_DEPTH);
        failed = !std::isfinite(candidate.association);
    } catch (...) {
        AD::abort_recording();
        failed = true;
    }
    bool recovered = false;
    try {
        const Candidate candidate = evaluate_candidate(valid_case, QUALIFYING_DEPTH);
        recovered = std::isfinite(candidate.association) && candidate.sites.size() == valid_case.parameters.site_count();
    } catch (...) {
        AD::abort_recording();
    }
    const bool passed = failed && recovered;
    recording_support::observe("failure-continuation", "status", "1", passed ? 1.0 : 0.0, std::nullopt,
        "invalid zero-total case followed by valid pure-2b case", "Diagnostic", "case_continuation",
        "failed case is contained and next case remains evaluable", "failure-continuation", "pure-2b", "status", "case", "failed_then_recovered", "case");
    return passed;
}

void association_comparison() {
    expect_keys();
    CppAD::ErrorHandler error_handler(cppad_error_handler);
    const auto cases = test_cases();
    const std::string selector = campaign_selector();
    const auto depths = campaign_depths(selector);
    std::vector<const TestCase*> selected;
    for (const auto& test_case : cases) {
        print_input_manifest(test_case);
        if (selected_case(test_case, selector)) selected.push_back(&test_case);
    }
    CHECK(!selected.empty());
    std::size_t sites = 0, gradient_entries = 0, hessian_entries = 0;
    std::size_t analytic_sites = 0;
    for (const TestCase* test_case : selected) {
        CHECK(test_case->amounts.size() == test_case->parameters.component_count());
        CHECK(std::all_of(test_case->amounts.begin(), test_case->amounts.end(), [](double value) { return std::isfinite(value) && value >= 0.0; }));
        double total = 0.0;
        for (double amount : test_case->amounts) total += amount;
        CHECK(total > 0.0);
        const std::size_t width = 2 + test_case->parameters.component_count() + (test_case->active.has_value() ? 1 : 0);
        sites += test_case->parameters.site_count();
        gradient_entries += width;
        hessian_entries += width * width;
        if (test_case->analytic) analytic_sites += test_case->parameters.site_count();
    }
    const auto expect_depth_group = [&](const char* group, const char* definition,
        std::size_t count, const char* coordinates) {
        for (const std::size_t depth : depths) {
            const std::string prefix = depth_prefix(depth, group);
            expect_group(recording_support::active_test_id(), group,
                prefix.c_str(), count, definition, coordinates);
        }
    };
    expect_depth_group("site-fraction", "physical site fractions", sites, "case,site");
    expect_depth_group("candidate-residual", "candidate residual", sites, "case,site");
    expect_depth_group("reference-residual", "reference residual", sites, "case,site");
    expect_depth_group("association-energy", "extensive association energy", selected.size(), "case");
    expect_depth_group("residual-energy", "extensive residual energy", selected.size(), "case");
    expect_depth_group("association-gradient", "scaled association gradient", gradient_entries, "case,coordinate");
    expect_depth_group("residual-gradient", "scaled residual gradient", gradient_entries, "case,coordinate");
    expect_depth_group("association-hessian", "scaled association Hessian", hessian_entries, "case,row,column");
    expect_depth_group("residual-hessian", "scaled residual Hessian", hessian_entries, "case,row,column");
    expect_depth_group("k-residual-identity", "shared K and residual identity", sites, "case,site");
    expect_depth_group("analytic-control", "pure 2B and 3B controls", analytic_sites, "case,site");
    expect_group(recording_support::active_test_id(), "reference-quality", "reference-quality/", selected.size(), "reference condition diagnostics", "case");
    expect_group(recording_support::active_test_id(), "reference-first-sensitivity", "reference-first-sensitivity/", selected.size(), "reference first sensitivity diagnostics", "coordinate");
    expect_group(recording_support::active_test_id(), "reference-second-sensitivity", "reference-second-sensitivity/", selected.size(), "reference second sensitivity diagnostics", "row,column");
    expect_group(recording_support::active_test_id(), "reference-sensitivity-proxy", "reference-sensitivity-proxy/", selected.size(), "reference sensitivity error proxy", "case");
    expect_group(recording_support::active_test_id(), "reference-newton-defect", "reference-newton-defect/", selected.size(), "reference Newton defect", "case");
    const bool full_campaign = selector == "all";
    if (full_campaign) {
        constexpr std::size_t centered_coordinates = 31;
        constexpr std::size_t centered_mixed = 10;
        expect_group(recording_support::active_test_id(), "centered-gradient", "centered-gradient/", 8 * centered_coordinates, "candidate and reference finite-difference gradients at two steps", "source,energy,case,step,coordinate");
        expect_group(recording_support::active_test_id(), "centered-hessian", "centered-hessian/", 8 * centered_coordinates, "candidate and reference finite-difference Hessians at two steps", "source,energy,case,step,coordinate");
        expect_group(recording_support::active_test_id(), "centered-mixed", "centered-mixed/", 8 * centered_mixed, "candidate and reference finite-difference mixed Hessians at two steps", "source,energy,case,step,row,column");
        expect_group(recording_support::active_test_id(), "production-weighted-hessian", "production-weighted-hessian/", 1, "known exact-zero production limitation", "case");
        expect_group(recording_support::active_test_id(), "extensive-scaling", "extensive-scaling/", 1, "extensive Helmholtz scaling", "case");
        expect_group(recording_support::active_test_id(), "failure-continuation", "failure-continuation", 1, "failed case followed by valid case", "case");
    }
    bool passed = true;
    for (const TestCase* test_case : selected) {
        Reference reference;
        bool reference_ok = true;
        try { reference = evaluate_reference(*test_case); }
        catch (const std::exception& error) {
            AD::abort_recording();
            reference_ok = false;
            reference.status = "reference_exception";
            reference.value_status = "reference_exception";
            reference.derivative_status = "reference_exception";
            reference.error = error.what();
        }
        const auto reference_inputs = inputs_for(*test_case, Candidate{}, reference, QUALIFYING_DEPTH);
        reference_diagnostics(*test_case, reference, reference_inputs);
        for (const std::size_t depth : depths) {
            state_status(*test_case, depth, "running");
            Candidate candidate;
            bool candidate_ok = true;
            try { candidate = evaluate_candidate(*test_case, depth); }
            catch (...) { AD::abort_recording(); candidate_ok = false; }
            const bool finite_candidate = candidate_ok && candidate.sites.size() == test_case->parameters.site_count()
                && std::isfinite(candidate.association) && std::isfinite(candidate.residual)
                && std::all_of(candidate.sites.begin(), candidate.sites.end(), [](double value) { return std::isfinite(value) && value > 0.0 && value <= 1.0; });
            bool state_passed = true;
            if (finite_candidate && reference_ok && reference.value_available) {
                state_passed = compare_values(*test_case, candidate, reference, depth) && state_passed;
                if (reference.derivative_available) {
                    state_passed = compare_derivative_set(*test_case, candidate, reference, candidate.association_gradient, candidate.association_hessian,
                        reference.association_tensor.gradient, reference.association_tensor.hessian, "association", "scaled association derivatives", depth) && state_passed;
                    state_passed = compare_derivative_set(*test_case, candidate, reference, candidate.residual_gradient, candidate.residual_hessian,
                        reference.tensor.gradient, reference.tensor.hessian, "residual", "scaled residual derivatives", depth) && state_passed;
                } else {
                    emit_derivative_unavailable(*test_case, candidate, reference, depth, reference.derivative_status.c_str());
                    state_passed = false;
                }
                state_passed = k_identity(*test_case, candidate, depth) && state_passed;
            } else {
                if (!reference_ok) reference.status = "reference_exception";
                if (!candidate_ok) reference.status = "candidate_exception";
                emit_unavailable(*test_case, candidate, reference, depth);
                state_passed = false;
            }
            state_passed = analytic_controls(*test_case, candidate, candidate_ok, depth) && state_passed;
            if (state_passed) state_status(*test_case, depth, "passed");
            else state_status(*test_case, depth, "criterion_or_reference_unavailable", reference.status.c_str());
            passed = state_passed && passed;
        }
        if (full_campaign) passed = production_zero_support_status(*test_case) && passed;
    }
    if (full_campaign) {
        passed = centered_checks(cases[0], {0, 1, 2}, {{0, 1}}) && passed;
        passed = centered_checks(cases[1], {0, 1, 2}, {{0, 1}}) && passed;
        passed = centered_checks(cases[3], {0, 1, 2}, {{0, 1}}) && passed;
        passed = centered_checks(cases[7], {0, 1, 2, 3}, {{0, 1}}) && passed;
        passed = centered_checks(cases[8], {4}, {{0, 4}}, {1.0e-3, 3.0e-4}) && passed;
        passed = centered_checks(cases[9], {4}, {{0, 4}}, {1.0e-3, 3.0e-4}) && passed;
        passed = centered_checks(cases[5], {0, 1, 2, 3}, {{0, 1}}) && passed;
        passed = centered_checks(cases[13], {0, 1, 2, 3}, {{0, 1}}) && passed;
        passed = centered_checks(cases[26], {0, 1, 2, 3}, {{0, 1}}) && passed;
        passed = centered_checks(cases[34], {0, 1, 2, 3}, {{0, 1}}) && passed;
        passed = extensive_scaling(cases[7]) && passed;
        passed = failure_continuation(cases[0]) && passed;
    }
    CHECK(passed);
}

} // namespace

int main(int argc, char** argv) {
    return recording_support::run(argc, argv, {
        {{"native.association.comparison", "association-comparison"}, association_comparison},
        {{"native.association.unit-ceiling", "association-unit-ceiling"}, association_unit_ceiling}
    });
}
