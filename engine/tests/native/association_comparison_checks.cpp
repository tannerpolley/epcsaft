#include <cppad/example/cppad_eigen.hpp>

#include "contributions/association.hpp"
#include "contributions/helmholtz.hpp"
#include "properties/helmholtz_tape.hpp"
#include "properties_support.hpp"
#include "recording_support.hpp"

#include <Eigen/LU>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
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
using Parameters = epcsaft::Parameters;
constexpr std::size_t DEPTH = 16;
constexpr double VALUE_TOL = 1.0e-10;
constexpr double SITE_RTOL = 1.0e-8;
constexpr double ENERGY_RTOL = 1.0e-9;
constexpr double REFERENCE_TOL = 1.0e-12;
constexpr double CONDITION_LIMIT = 1.0e12;
constexpr double DERIVATIVE_ATOL = 1.0e-8;
constexpr double DERIVATIVE_RTOL = 1.0e-6;
constexpr double CENTERED_ATOL = 1.0e-7;
constexpr double CENTERED_RTOL = 1.0e-5;

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
struct Reference final {
    std::vector<double> sites;
    std::vector<double> mass_action;
    double association = std::numeric_limits<double>::quiet_NaN();
    double residual = std::numeric_limits<double>::quiet_NaN();
    properties_support::ReferenceTensor tensor;
    properties_support::ReferenceTensor association_tensor;
    double condition = std::numeric_limits<double>::quiet_NaN();
    double linear_residual = std::numeric_limits<double>::quiet_NaN();
    std::string status = "unavailable";
    bool available = false;
};

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
    const ParameterSet& parameters, const std::vector<Scalar>& delta) {
    const std::size_t sites = parameters.site_component.size();
    const auto k = association_k(density, composition, parameters, delta);
    std::vector<Scalar> y(sites, Scalar(0.0));
    for (std::size_t row = 0; row < sites; ++row) {
        Scalar row_sum = Scalar(0.0);
        for (std::size_t column = 0; column < sites; ++column) row_sum += k[row * sites + column];
        y[row] = epcsaft::scalar_log(Scalar(2.0) / (Scalar(1.0)
            + epcsaft::scalar_sqrt(Scalar(1.0) + Scalar(4.0) * row_sum)));
    }
    for (std::size_t iteration = 0; iteration < DEPTH; ++iteration) {
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
Tape make_tape(const TestCase& test_case) {
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
    const auto sites = fixed_candidate(density, mole_fractions, parameters, delta);
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
Candidate evaluate_candidate(const TestCase& test_case) {
    Tape tape = make_tape(test_case);
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
Reference evaluate_reference(const TestCase& test_case) {
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
    const auto k = association_k(density, mole_fractions, test_case.parameters, delta);
    const std::size_t sites = result.sites.size();
    Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(sites), static_cast<Eigen::Index>(sites));
    for (std::size_t row = 0; row < sites; ++row) {
        double bonded = 0.0;
        for (std::size_t column = 0; column < sites; ++column) bonded += k[row * sites + column] * result.sites[column];
        for (std::size_t column = 0; column < sites; ++column)
            jacobian(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) =
                (row == column ? 1.0 + bonded : 0.0) + result.sites[row] * k[row * sites + column];
    }
    bool quality = std::all_of(result.mass_action.begin(), result.mass_action.end(),
        [](double value) { return std::isfinite(value); });
    double residual_norm = 0.0;
    for (double value : result.mass_action) residual_norm = std::max(residual_norm, std::abs(value));
    if (!quality || residual_norm > REFERENCE_TOL) { quality = false; result.status = "reference_residual_failure"; }
    if (quality && sites != 0) {
        Eigen::FullPivLU<Eigen::MatrixXd> factorization(jacobian);
        if (factorization.rank() != jacobian.rows()) {
            quality = false;
            result.status = "reference_rank_failure";
        } else {
            const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(jacobian.rows(), jacobian.cols());
            const Eigen::MatrixXd inverse = factorization.solve(identity);
            result.linear_residual = (jacobian * inverse - identity).cwiseAbs().maxCoeff();
            result.condition = jacobian.norm() * inverse.norm();
            if (!std::isfinite(result.linear_residual) || !std::isfinite(result.condition)
                || result.linear_residual > 1.0e-10 || result.condition > CONDITION_LIMIT) {
                quality = false;
                result.status = "reference_condition_failure";
            }
        }
    } else if (quality) {
        result.condition = 1.0;
        result.linear_residual = 0.0;
    }
    if (quality) {
        try {
            result.tensor = properties_support::reference_tensor(test_case.parameters, test_case.temperature,
                test_case.volume, test_case.amounts, result.sites, test_case.active);
            result.association_tensor = properties_support::reference_tensor(test_case.parameters,
                test_case.temperature, test_case.volume, test_case.amounts, result.sites, test_case.active, true);
            result.available = std::all_of(result.tensor.gradient.begin(), result.tensor.gradient.end(),
                [](double value) { return std::isfinite(value); })
                && std::all_of(result.tensor.hessian.begin(), result.tensor.hessian.end(),
                    [](double value) { return std::isfinite(value); })
                && std::all_of(result.association_tensor.gradient.begin(), result.association_tensor.gradient.end(),
                    [](double value) { return std::isfinite(value); })
                && std::all_of(result.association_tensor.hessian.begin(), result.association_tensor.hessian.end(),
                    [](double value) { return std::isfinite(value); });
            if (!result.available) { quality = false; result.status = "reference_tensor_nonfinite"; }
        } catch (const std::exception&) {
            AD::abort_recording();
            quality = false;
            result.status = "reference_tensor_failure";
        }
    }
    if (quality) result.status = "passed";
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
    bool analytic = false, std::optional<std::uint32_t> seed = std::nullopt) {
    return {id, std::move(parameters), temperature, volume, std::move(amounts), active, analytic, seed};
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

std::vector<TestCase> test_cases() {
    auto weak = pure_two_b();
    weak.association_energy_over_k[1] = weak.association_energy_over_k[2] = 350.0;
    weak.association_volume[1] = weak.association_volume[2] = 0.01;
    auto strong = pure_two_b();
    strong.association_energy_over_k[1] = strong.association_energy_over_k[2] = 3200.0;
    strong.association_volume[1] = strong.association_volume[2] = 0.20;
    const auto associating = properties_support::associating_binary();
    return {
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
}

std::string active_name(const TestCase& test_case) {
    if (!test_case.active.has_value()) return "none";
    return test_case.active->family == epcsaft::ParameterFamily::AssociationEnergyOverK
        ? "association_energy_over_k" : "association_volume";
}

std::string inputs_for(const TestCase& test_case, const Candidate& candidate, const Reference& reference) {
    return recording_support::fields({recording_support::field("case_id", test_case.id),
        recording_support::field("T_K", test_case.temperature), recording_support::field("V_m3", test_case.volume),
        recording_support::field("amounts_mol", test_case.amounts), recording_support::field("active_parameter", active_name(test_case)),
        recording_support::field("generation_seed", test_case.seed.has_value() ? std::to_string(*test_case.seed) : "none"),
        recording_support::field("candidate_sites", candidate.sites), recording_support::field("reference_sites", reference.sites),
        recording_support::field("reference_status", reference.status)});
}

double scale_for(const TestCase& test_case, std::size_t coordinate) {
    const std::size_t components = test_case.parameters.component_count();
    if (coordinate == 0) return 300.0;
    if (coordinate == 1) return 1.0e-4;
    if (coordinate < 2 + components) return 1.0;
    return test_case.active->family == epcsaft::ParameterFamily::AssociationEnergyOverK ? 1000.0 : 0.1;
}

void reference_diagnostics(const TestCase& test_case, const Reference& reference, const std::string& inputs) {
    recording_support::observe((std::string("reference-quality/") + test_case.id).c_str(), "reference_condition", "1",
        reference.condition, reference.condition, inputs, "Diagnostic", "implicit_reference_quality", reference.status.c_str(),
        "reference-quality", test_case.id, "condition", "T,V,n,active", "J condition", "case");
    recording_support::observe((std::string("reference-linear-residual/") + test_case.id).c_str(), "reference_linear_residual", "1",
        reference.linear_residual, reference.linear_residual, inputs, "Diagnostic", "implicit_reference_quality", reference.status.c_str(),
        "reference-quality", test_case.id, "linear_residual", "T,V,n,active", "J J^-1-I", "case");
}

bool compare_row(const char* comparison_id, const char* quantity, const char* unit, double model, double reference,
    const std::string& inputs, const char* basis, const char* method, const char* claim, const char* criterion,
    double atol, double rtol, bool, const char* comparison_group_id = nullptr, const char* point_id = nullptr,
    const char* quantity_latex = nullptr, const char* inputs_latex = nullptr, const char* relationship_latex = nullptr,
    const char* coordinates = nullptr) {
    const bool verdict = recording_support::scalar_verdict(model, reference, criterion, atol, rtol);
    recording_support::print_json_compare(comparison_id, quantity, unit, model, reference, inputs, basis, method, claim,
        criterion, atol, rtol, verdict, verdict ? "Criterion satisfied" : "Criterion failed", comparison_group_id,
        point_id, quantity_latex, inputs_latex, relationship_latex, coordinates);
    return verdict;
}

bool compare_values(const TestCase& test_case, const Candidate& candidate, const Reference& reference) {
    bool passed = true;
    const auto inputs = inputs_for(test_case, candidate, reference);
    reference_diagnostics(test_case, reference, inputs);
    for (std::size_t site = 0; site < candidate.sites.size(); ++site) {
        const std::string point = std::string(test_case.id) + "/" + std::to_string(site);
        passed = compare_row((std::string("site-fraction/") + point).c_str(), "site_fraction", "1",
            candidate.sites[site], reference.sites[site], inputs, "Comparison", "fixed_depth_candidate", "site fraction agreement",
            "reference_scale", 1.0e-12, SITE_RTOL, true, "site-fraction", point.c_str(), "X", "T,V,n", "X_candidate=X_reference", "site") && passed;
        passed = compare_row((std::string("candidate-residual/") + point).c_str(), "mass_action_residual", "1",
            candidate.mass_action[site], 0.0, inputs, "Comparison", "defining_residual", "candidate original residual",
            "strict_absolute", VALUE_TOL, 0.0, true, "candidate-residual", point.c_str(), "F", "T,V,n", "F_candidate=0", "site") && passed;
        passed = compare_row((std::string("reference-residual/") + point).c_str(), "reference_mass_action_residual", "1",
            reference.mass_action[site], 0.0, inputs, "Reference", "defining_residual", "reference original residual",
            "strict_absolute", REFERENCE_TOL, 0.0, true, "reference-residual", point.c_str(), "F", "T,V,n", "F_reference=0", "site") && passed;
    }
    passed = compare_row((std::string("association-energy/") + test_case.id).c_str(), "Phi_association", "1",
        candidate.association, reference.association, inputs, "Comparison", "fixed_depth_candidate", "extensive association energy",
        "strict_reference_or_one", 0.0, ENERGY_RTOL, true, "association-energy", test_case.id, "Phi_assoc", "T,V,n", "candidate=reference", "case") && passed;
    passed = compare_row((std::string("residual-energy/") + test_case.id).c_str(), "Phi_residual", "1",
        candidate.residual, reference.residual, inputs, "Comparison", "fixed_depth_candidate", "extensive residual energy",
        "strict_reference_or_one", 0.0, ENERGY_RTOL, true, "residual-energy", test_case.id, "Phi_res", "T,V,n", "candidate=reference", "case") && passed;
    return passed;
}

bool compare_derivative_set(const TestCase& test_case, const Candidate& candidate, const Reference& reference,
    const std::vector<double>& model_gradient, const std::vector<double>& model_hessian,
    const std::vector<double>& reference_gradient, const std::vector<double>& reference_hessian,
    const char* label, const char* claim) {
    const std::size_t width = model_gradient.size();
    const auto inputs = inputs_for(test_case, candidate, reference);
    bool passed = true;
    const std::string gradient_group = std::string(label) + "-gradient";
    const std::string hessian_group = std::string(label) + "-hessian";
    for (std::size_t coordinate = 0; coordinate < width; ++coordinate) {
        const double scale = scale_for(test_case, coordinate);
        const std::string point = std::string(test_case.id) + "/" + std::to_string(coordinate);
        const auto details = recording_support::fields({inputs, recording_support::field("coordinate", coordinate),
            recording_support::field("scale", scale), recording_support::field("model_raw", model_gradient[coordinate]),
            recording_support::field("reference_raw", reference_gradient[coordinate])});
        passed = compare_row((gradient_group + "/" + point).c_str(), "scaled_gradient", "1",
            model_gradient[coordinate] * scale, reference_gradient[coordinate] * scale, details, "Comparison", "CppAD_vs_implicit_reference",
            claim, "reference_scale", DERIVATIVE_ATOL, DERIVATIVE_RTOL, true, gradient_group.c_str(), point.c_str(), "dPhi/dq_scaled", "T,V,n,active", "candidate=reference", "coordinate") && passed;
    }
    for (std::size_t row = 0; row < width; ++row) for (std::size_t column = 0; column < width; ++column) {
        const std::size_t entry = row * width + column;
        const double scale = scale_for(test_case, row) * scale_for(test_case, column);
        const std::string point = std::string(test_case.id) + "/" + std::to_string(row) + "-" + std::to_string(column);
        const auto details = recording_support::fields({inputs, recording_support::field("row", row), recording_support::field("column", column),
            recording_support::field("scale", scale), recording_support::field("model_raw", model_hessian[entry]),
            recording_support::field("reference_raw", reference_hessian[entry])});
            passed = compare_row((hessian_group + "/" + point).c_str(), "scaled_hessian", "1",
            model_hessian[entry] * scale, reference_hessian[entry] * scale, details, "Comparison", "CppAD_vs_implicit_reference",
            claim, "reference_scale", DERIVATIVE_ATOL, DERIVATIVE_RTOL, true, hessian_group.c_str(), point.c_str(), "d2Phi/dq_scaled2", "T,V,n,active", "candidate=reference", "row,column") && passed;
    }
    return passed;
}

double candidate_value(const TestCase& test_case) {
    try {
        Tape tape = make_tape(test_case);
        return tape.function.Forward(0, tape.point)[1];
    } catch (...) {
        AD::abort_recording();
        return std::numeric_limits<double>::quiet_NaN();
    }
}

void emit_centered_unavailable(const TestCase& test_case, double relative_step, bool mixed) {
    const std::size_t width = 2 + test_case.parameters.component_count();
    for (std::size_t coordinate = 0; coordinate < width; ++coordinate) {
        const std::string point = std::string(test_case.id) + "/h" + recording_support::number(relative_step) + "/" + std::to_string(coordinate);
        const auto inputs = recording_support::fields({recording_support::field("case_id", test_case.id), recording_support::field("coordinate", coordinate), recording_support::field("status", "center_unavailable")});
        recording_support::observe((std::string("centered-gradient/") + point).c_str(), "centered_gradient", "1", std::numeric_limits<double>::quiet_NaN(), std::nullopt, inputs, "Independent_check", "five_point_candidate", "candidate center gradient", "centered-gradient", point.c_str(), "dPhi/dq_scaled", "T,V,n", "status", "coordinate");
        recording_support::observe((std::string("centered-hessian/") + point).c_str(), "centered_hessian", "1", std::numeric_limits<double>::quiet_NaN(), std::nullopt, inputs, "Independent_check", "five_point_candidate", "candidate center Hessian", "centered-hessian", point.c_str(), "d2Phi/dq_scaled2", "T,V,n", "status", "coordinate");
    }
    if (mixed && width > 1) recording_support::observe((std::string("centered-mixed/") + test_case.id + "/0-1").c_str(), "centered_mixed_hessian", "1", std::numeric_limits<double>::quiet_NaN(), std::nullopt, recording_support::fields({recording_support::field("case_id", test_case.id), recording_support::field("status", "center_unavailable")}), "Independent_check", "four_point_candidate", "candidate center mixed Hessian", "centered-mixed", test_case.id, "d2Phi/dTdV", "T,V,n", "status", "row,column");
}

bool centered_checks(const TestCase& test_case) {
    bool passed = true;
    for (const double relative_step : {1.0e-4, 3.0e-5}) {
        Tape tape;
        std::vector<double> center, jacobian, hessian;
        try {
            tape = make_tape(test_case);
            center = tape.function.Forward(0, tape.point);
            jacobian = tape.function.Jacobian(tape.point);
            hessian = tape.function.Hessian(tape.point, 1);
        } catch (...) {
            AD::abort_recording();
            emit_centered_unavailable(test_case, relative_step, relative_step == 1.0e-4);
            passed = false;
            continue;
        }
        for (std::size_t coordinate = 0; coordinate < tape.width; ++coordinate) {
            const double scale = scale_for(test_case, coordinate);
            const double step = relative_step * scale;
            auto shifted = [&](double amount) {
                TestCase result = test_case;
                if (coordinate == 0) result.temperature += amount;
                else if (coordinate == 1) result.volume += amount;
                else result.amounts[coordinate - 2] += amount;
                return result;
            };
            const double m2 = candidate_value(shifted(-2.0 * step));
            const double m1 = candidate_value(shifted(-step));
            const double p1 = candidate_value(shifted(step));
            const double p2 = candidate_value(shifted(2.0 * step));
            const double first = (m2 - 8.0 * m1 + 8.0 * p1 - p2) / (12.0 * step);
            const double second = (-m2 + 16.0 * m1 - 30.0 * center[1] + 16.0 * p1 - p2) / (12.0 * step * step);
            const std::string point = std::string(test_case.id) + "/h" + recording_support::number(relative_step) + "/" + std::to_string(coordinate);
            const auto inputs = recording_support::fields({recording_support::field("case_id", test_case.id),
                recording_support::field("coordinate", coordinate), recording_support::field("scale", scale), recording_support::field("step", step),
                recording_support::field("center_candidate", center[1])});
            passed = compare_row((std::string("centered-gradient/") + point).c_str(), "centered_gradient", "1",
                first * scale, jacobian[tape.width + coordinate] * scale, inputs, "Independent_check", "five_point_candidate", "candidate center gradient",
                "reference_scale", CENTERED_ATOL, CENTERED_RTOL, true, "centered-gradient", point.c_str(), "dPhi/dq_scaled", "T,V,n", "stencil=graph", "coordinate") && passed;
            passed = compare_row((std::string("centered-hessian/") + point).c_str(), "centered_hessian", "1",
                second * scale * scale, hessian[coordinate * tape.width + coordinate] * scale * scale, inputs, "Independent_check", "five_point_candidate", "candidate center Hessian",
                "reference_scale", CENTERED_ATOL, CENTERED_RTOL, true, "centered-hessian", point.c_str(), "d2Phi/dq_scaled2", "T,V,n", "stencil=graph", "coordinate") && passed;
        }
        if (relative_step == 1.0e-4 && tape.width > 1) {
            const double t_step = relative_step * scale_for(test_case, 0);
            const double v_step = relative_step * scale_for(test_case, 1);
            const auto shift = [&](double dt, double dv) {
                TestCase result = test_case;
                result.temperature += dt;
                result.volume += dv;
                return result;
            };
            const double mixed = (candidate_value(shift(t_step, v_step)) - candidate_value(shift(t_step, -v_step))
                - candidate_value(shift(-t_step, v_step)) + candidate_value(shift(-t_step, -v_step))) / (4.0 * t_step * v_step);
            passed = compare_row((std::string("centered-mixed/") + test_case.id + "/0-1").c_str(), "centered_mixed_hessian", "1",
                mixed * scale_for(test_case, 0) * scale_for(test_case, 1),
                tape.function.Hessian(tape.point, 1)[1] * scale_for(test_case, 0) * scale_for(test_case, 1),
                recording_support::fields({recording_support::field("case_id", test_case.id), recording_support::field("row", 0), recording_support::field("column", 1),
                    recording_support::field("T_step", t_step), recording_support::field("V_step", v_step), recording_support::field("center_candidate", center[1])}),
                "Independent_check", "four_point_candidate", "candidate center mixed Hessian", "reference_scale", CENTERED_ATOL, CENTERED_RTOL, true,
                "centered-mixed", test_case.id, "d2Phi/dTdV", "T,V,n", "stencil=graph", "row,column") && passed;
        }
    }
    return passed;
}

bool analytic_controls(const TestCase& test_case, const Candidate& candidate, bool candidate_ok) {
    if (!test_case.analytic) return true;
    const std::size_t expected_size = test_case.parameters.site_count();
    if (!candidate_ok || candidate.sites.size() != expected_size) {
        for (std::size_t site = 0; site < expected_size; ++site)
            recording_support::observe((std::string("analytic-control/") + test_case.id + "/" + std::to_string(site)).c_str(),
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
        passed = compare_row((std::string("analytic-control/") + test_case.id + "/" + std::to_string(site)).c_str(),
            "site_fraction", "1", candidate.sites[site], expected[site], recording_support::fields({recording_support::field("case_id", test_case.id), recording_support::field("coupling", k[1])}),
            "Analytic_control", "closed_form", "candidate matches pure-site control", "strict_reference_scale", 1.0e-12, 1.0e-12, true,
            "analytic-control", test_case.id, "X", "T,V,n", "candidate=closed_form", "site") && passed;
    }
    return passed;
}

bool k_identity(const TestCase& test_case, const Candidate& candidate) {
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
        passed = compare_row((std::string("k-residual-identity/") + test_case.id + "/" + std::to_string(row)).c_str(),
            "mass_action_residual", "1", candidate.mass_action[row], manual,
            recording_support::fields({recording_support::field("case_id", test_case.id), recording_support::field("K_row", std::vector<double>(k.begin() + row * candidate.sites.size(), k.begin() + (row + 1) * candidate.sites.size()))}),
            "Relationship", "shared_residual", "association residual uses rho_hat x_b nu_b Delta", "strict_reference_scale", 1.0e-14, 1.0e-14, true,
            "k-residual-identity", test_case.id, "F", "T,V,n", "F_shared=F_manual", "site") && passed;
    }
    return passed;
}

void emit_unavailable(const TestCase& test_case, const Candidate& candidate, const Reference& reference) {
    const auto inputs = inputs_for(test_case, candidate, reference);
    const std::size_t sites = test_case.parameters.site_count();
    const std::size_t width = 2 + test_case.parameters.component_count() + (test_case.active.has_value() ? 1 : 0);
    reference_diagnostics(test_case, reference, inputs);
    for (std::size_t site = 0; site < sites; ++site) {
        const auto candidate_site = site < candidate.sites.size() ? std::optional<double>(candidate.sites[site]) : std::nullopt;
        const auto reference_site = site < reference.sites.size() ? std::optional<double>(reference.sites[site]) : std::nullopt;
        recording_support::observe((std::string("site-fraction/") + test_case.id + "/" + std::to_string(site)).c_str(), "site_fraction", "1",
            candidate_site.value_or(std::numeric_limits<double>::quiet_NaN()), reference_site, inputs, "Comparison", "fixed_depth_candidate", reference.status.c_str(), "site-fraction", test_case.id, "X", "T,V,n", "status", "site");
        recording_support::observe((std::string("candidate-residual/") + test_case.id + "/" + std::to_string(site)).c_str(), "mass_action_residual", "1",
            site < candidate.mass_action.size() ? candidate.mass_action[site] : std::numeric_limits<double>::quiet_NaN(), std::nullopt, inputs, "Comparison", "defining_residual", reference.status.c_str(), "candidate-residual", test_case.id, "F", "T,V,n", "status", "site");
        recording_support::observe((std::string("reference-residual/") + test_case.id + "/" + std::to_string(site)).c_str(), "reference_mass_action_residual", "1",
            site < reference.mass_action.size() ? reference.mass_action[site] : std::numeric_limits<double>::quiet_NaN(), std::nullopt, inputs, "Reference", "defining_residual", reference.status.c_str(), "reference-residual", test_case.id, "F", "T,V,n", "status", "site");
    }
    recording_support::observe((std::string("association-energy/") + test_case.id).c_str(), "Phi_association", "1", candidate.association, std::nullopt, inputs, "Comparison", "fixed_depth_candidate", reference.status.c_str(), "association-energy", test_case.id, "Phi_assoc", "T,V,n", "status", "case");
    recording_support::observe((std::string("residual-energy/") + test_case.id).c_str(), "Phi_residual", "1", candidate.residual, std::nullopt, inputs, "Comparison", "fixed_depth_candidate", reference.status.c_str(), "residual-energy", test_case.id, "Phi_res", "T,V,n", "status", "case");
    for (std::size_t coordinate = 0; coordinate < width; ++coordinate) for (const char* label : {"association", "residual"})
        recording_support::observe((std::string(label) + "-gradient/" + test_case.id + "/" + std::to_string(coordinate)).c_str(), "scaled_gradient", "1", std::numeric_limits<double>::quiet_NaN(), std::nullopt, inputs, "Comparison", "CppAD_vs_implicit_reference", reference.status.c_str(), (std::string(label) + "-gradient").c_str(), test_case.id, "dPhi/dq", "T,V,n,active", "status", "coordinate");
    for (std::size_t entry = 0; entry < width * width; ++entry) for (const char* label : {"association", "residual"})
        recording_support::observe((std::string(label) + "-hessian/" + test_case.id + "/" + std::to_string(entry)).c_str(), "scaled_hessian", "1", std::numeric_limits<double>::quiet_NaN(), std::nullopt, inputs, "Comparison", "CppAD_vs_implicit_reference", reference.status.c_str(), (std::string(label) + "-hessian").c_str(), test_case.id, "d2Phi/dq2", "T,V,n,active", "status", "entry");
    for (std::size_t site = 0; site < sites; ++site)
        recording_support::observe((std::string("k-residual-identity/") + test_case.id + "/" + std::to_string(site)).c_str(),
            "mass_action_residual", "1", std::numeric_limits<double>::quiet_NaN(), std::nullopt, inputs,
            "Relationship", "shared_residual", reference.status.c_str(), "k-residual-identity", test_case.id, "F", "T,V,n", "status", "site");
}

bool production_zero_support_status(const TestCase& test_case) {
    if (std::string(test_case.id) != "exact-zero-support") return true;
    bool failed_as_expected = false;
    std::string error_message;
    try {
        epcsaft::HelmholtzModel model(test_case.parameters);
        model.evaluate(test_case.temperature, test_case.volume, test_case.amounts, epcsaft::DerivativeOrder::Hessian);
    } catch (const std::exception& error) {
        error_message = error.what();
        failed_as_expected = error_message.find("association Hessian is singular") != std::string::npos;
    }
    recording_support::observe("production-weighted-hessian/exact-zero-support", "status", "1",
        failed_as_expected ? 1.0 : 0.0, std::nullopt,
        "case_id:exact-zero-support;known_status:association Hessian is singular;exception=" + error_message,
        "Diagnostic", "production_probe", "known exact-zero support limitation", "production-weighted-hessian", "exact-zero-support", "status", "T,V,n", "failure expected", "case");
    return failed_as_expected;
}

bool extensive_scaling(const TestCase& test_case) {
    TestCase scaled = test_case;
    for (double& amount : scaled.amounts) amount *= 2.0;
    scaled.volume *= 2.0;
    const double center = candidate_value(test_case);
    const double doubled = candidate_value(scaled);
    return compare_row("extensive-scaling/heterogeneous-cross", "Phi_residual", "1", doubled, 2.0 * center,
        "case_id:heterogeneous-cross;amount_scale:2;volume_scale:2", "Relationship", "candidate_extensivity", "Phi_res scales with total amount at fixed density",
        "strict_reference_scale", 1.0e-10, 1.0e-12, true, "extensive-scaling", "heterogeneous-cross", "Phi_res", "T,V,n", "Phi(2n,2V)=2Phi(n,V)", "case");
}

bool failure_continuation(const TestCase& valid_case) {
    TestCase invalid = valid_case;
    std::fill(invalid.amounts.begin(), invalid.amounts.end(), 0.0);
    bool failed = false;
    try {
        const Candidate candidate = evaluate_candidate(invalid);
        failed = !std::isfinite(candidate.association);
    } catch (...) {
        AD::abort_recording();
        failed = true;
    }
    bool recovered = false;
    try {
        const Candidate candidate = evaluate_candidate(valid_case);
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
    CppAD::ErrorHandler error_handler(cppad_error_handler);
    const auto cases = test_cases();
    std::size_t sites = 0, gradient_entries = 0, hessian_entries = 0;
    for (const auto& test_case : cases) {
        CHECK(test_case.amounts.size() == test_case.parameters.component_count());
        CHECK(std::all_of(test_case.amounts.begin(), test_case.amounts.end(), [](double value) { return std::isfinite(value) && value >= 0.0; }));
        double total = 0.0;
        for (double amount : test_case.amounts) total += amount;
        CHECK(total > 0.0);
        const std::size_t width = 2 + test_case.parameters.component_count() + (test_case.active.has_value() ? 1 : 0);
        sites += test_case.parameters.site_count();
        gradient_entries += width;
        hessian_entries += width * width;
    }
    const auto& centered_pure = cases[0];
    const auto& centered = cases[7];
    const std::size_t centered_width = 2 + centered.parameters.component_count();
    const std::size_t centered_pure_width = 2 + centered_pure.parameters.component_count();
    recording_support::expect_group(recording_support::active_test_id(), "site-fraction", "site-fraction/", sites, "physical site fractions", "case,site");
    recording_support::expect_group(recording_support::active_test_id(), "candidate-residual", "candidate-residual/", sites, "candidate residual", "case,site");
    recording_support::expect_group(recording_support::active_test_id(), "reference-residual", "reference-residual/", sites, "reference residual", "case,site");
    recording_support::expect_group(recording_support::active_test_id(), "association-energy", "association-energy/", cases.size(), "extensive association energy", "case");
    recording_support::expect_group(recording_support::active_test_id(), "residual-energy", "residual-energy/", cases.size(), "extensive residual energy", "case");
    recording_support::expect_group(recording_support::active_test_id(), "association-gradient", "association-gradient/", gradient_entries, "scaled association gradient", "case,coordinate");
    recording_support::expect_group(recording_support::active_test_id(), "residual-gradient", "residual-gradient/", gradient_entries, "scaled residual gradient", "case,coordinate");
    recording_support::expect_group(recording_support::active_test_id(), "association-hessian", "association-hessian/", hessian_entries, "scaled association Hessian", "case,row,column");
    recording_support::expect_group(recording_support::active_test_id(), "residual-hessian", "residual-hessian/", hessian_entries, "scaled residual Hessian", "case,row,column");
    recording_support::expect_group(recording_support::active_test_id(), "reference-quality", "reference-quality/", cases.size(), "reference condition diagnostics", "case");
    recording_support::expect_group(recording_support::active_test_id(), "reference-linear-residual", "reference-linear-residual/", cases.size(), "reference linear residual diagnostics", "case");
    recording_support::expect_group(recording_support::active_test_id(), "centered-gradient", "centered-gradient/", 2 * (centered_width + centered_pure_width), "candidate centered gradients at two steps", "case,step,coordinate");
    recording_support::expect_group(recording_support::active_test_id(), "centered-hessian", "centered-hessian/", 2 * (centered_width + centered_pure_width), "candidate centered diagonal Hessians at two steps", "case,step,coordinate");
    recording_support::expect_group(recording_support::active_test_id(), "centered-mixed", "centered-mixed/", 2, "candidate centered mixed Hessian", "case,row,column");
    recording_support::expect_group(recording_support::active_test_id(), "analytic-control", "analytic-control/", 5, "pure 2B and 3B controls", "case,site");
    recording_support::expect_group(recording_support::active_test_id(), "k-residual-identity", "k-residual-identity/", sites, "shared K and residual identity", "case,site");
    recording_support::expect_group(recording_support::active_test_id(), "production-weighted-hessian", "production-weighted-hessian/", 1, "known exact-zero production limitation", "case");
    recording_support::expect_group(recording_support::active_test_id(), "extensive-scaling", "extensive-scaling/", 1, "extensive Helmholtz scaling", "case");
    recording_support::expect_group(recording_support::active_test_id(), "failure-continuation", "failure-continuation", 1, "failed case followed by valid case", "case");
    bool passed = true;
    for (const auto& test_case : cases) {
        Candidate candidate;
        Reference reference;
        bool candidate_ok = true;
        bool reference_ok = true;
        try { candidate = evaluate_candidate(test_case); }
        catch (...) { AD::abort_recording(); candidate_ok = false; }
        try { reference = evaluate_reference(test_case); }
        catch (...) { AD::abort_recording(); reference_ok = false; reference.status = "reference_exception"; }
        const bool finite_candidate = candidate_ok && candidate.sites.size() == test_case.parameters.site_count()
            && std::isfinite(candidate.association) && std::isfinite(candidate.residual)
            && std::all_of(candidate.sites.begin(), candidate.sites.end(), [](double value) { return std::isfinite(value) && value > 0.0 && value <= 1.0; });
        if (finite_candidate && reference_ok && reference.available) {
            passed = compare_values(test_case, candidate, reference) && passed;
            passed = compare_derivative_set(test_case, candidate, reference, candidate.association_gradient, candidate.association_hessian,
                reference.association_tensor.gradient, reference.association_tensor.hessian, "association", "scaled association derivatives") && passed;
            passed = compare_derivative_set(test_case, candidate, reference, candidate.residual_gradient, candidate.residual_hessian,
                reference.tensor.gradient, reference.tensor.hessian, "residual", "scaled residual derivatives") && passed;
            passed = k_identity(test_case, candidate) && passed;
        } else {
            if (!reference_ok) reference.status = "reference_exception";
            if (!candidate_ok) reference.status = "candidate_exception";
            emit_unavailable(test_case, candidate, reference);
            passed = false;
        }
        passed = analytic_controls(test_case, candidate, candidate_ok) && passed;
        passed = production_zero_support_status(test_case) && passed;
    }
    passed = centered_checks(centered_pure) && passed;
    passed = centered_checks(centered) && passed;
    passed = extensive_scaling(cases[7]) && passed;
    passed = failure_continuation(cases[0]) && passed;
    CHECK(passed);
}

} // namespace

int main(int argc, char** argv) {
    return recording_support::run(argc, argv, {
        {{"native.association.comparison", "Fixed_depth_native_association_comparison", "Numerical", "Contributions/Properties", "Explicit", "Numerical", "Association", "P2/M4", "site-fraction/*,candidate-residual/*,reference-residual/*,association-energy/*,residual-energy/*,association-gradient/*,residual-gradient/*,association-hessian/*,residual-hessian/*,reference-quality/*,reference-linear-residual/*,centered-gradient/*,centered-hessian/*,centered-mixed/*,analytic-control/*,k-residual-identity/*,production-weighted-hessian/*,extensive-scaling/*,failure-continuation", "engine/tests/native/association_comparison_checks.cpp", "association_comparison_checks", "association-comparison", ""}, association_comparison}
    });
}
