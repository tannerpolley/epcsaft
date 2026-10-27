#pragma once

#include "core/numerical_policy.hpp"
#include "core/trace.hpp"
#include <cstddef>
#include <cmath>
#include <limits>
#include <string>

namespace epcsaft::core {

template <typename Factor>
NumericalRank numerical_rank(const Factor& lu, const char* basis) {
    return {basis, static_cast<std::size_t>(lu.cols()), static_cast<std::size_t>(lu.rank()), double(lu.threshold())};
}

template <typename Factor>
void retain_rank(std::vector<NumericalRank>& records, const Factor& factor, const char* basis) {
    const auto rank = numerical_rank(factor, basis);
    retain_rank(records, rank);
    Scope scope("numerical_rank");
    if (scope.diagnostics_on()) scope.probe("rank", {double(rank.rank), double(rank.dimension), rank.threshold}, 0, basis);
}

// Rigal-Gaches normwise backward error; the zero right-hand side needs no unit floor.
inline double backward_error(double residual, double jacobian_norm, double solution, double right) {
    const double denominator = jacobian_norm * solution + right;
    return residual == 0.0 ? 0.0 : denominator > 0.0 ? residual / denominator : std::numeric_limits<double>::infinity();
}

// Ogita-Rump-Oishi Dot2: FMA TwoProduct and TwoSum preserve the cancellation defect.
template <typename Matrix, typename Solution, typename Right>
typename Right::PlainObject compensated_defect(const Matrix& matrix, const Solution& solution, const Right& right) {
    typename Right::PlainObject result(right.size());
    for (std::ptrdiff_t row = 0; row < right.size(); ++row) {
        double sum = right(row), error = 0.0;
        for (std::ptrdiff_t column = 0; column < solution.size(); ++column) {
            const double term = -matrix(row, column) * solution(column), next = sum + term, part = next - sum;
            error += (sum - (next - part)) + (term - part) + std::fma(-matrix(row, column), solution(column), -term);
            sum = next;
        }
        result(row) = sum + error;
    }
    return result;
}

// Undo numerical equilibration before measuring the supplied matrix's inverse norm.
template <typename Factor, typename Vector>
auto inverse_infinity_norm(const Factor& lu, const Vector& rows, const Vector& columns) {
    Vector sums = Vector::Zero(columns.size());
    for (std::ptrdiff_t i = 0; i < columns.size(); ++i) {
        const Vector rhs = Vector::Unit(columns.size(), i);
        sums += (columns.asDiagonal() * lu.solve(rows.asDiagonal() * rhs)).cwiseAbs();
    }
    return sums.maxCoeff();
}

// Propagate nonnegative row errors through an already admitted factorization.
template <typename Factor, typename Right>
typename Right::PlainObject inverse_error_action(const Factor& lu, const Right& right) {
    using Vector = typename Right::PlainObject;
    Vector result = Vector::Zero(right.size());
    for (std::ptrdiff_t i = 0; i < right.size(); ++i)
        result += lu.solve(Vector::Unit(right.size(), i)).cwiseAbs() * right[i];
    return result;
}

template <typename Real>
Real reduced_error_envelope(Real inverse, Real residual, Real floor, Real displacement = 0) {
    return inverse * (residual + floor) + displacement;
}

template <typename Matrix, typename Scale, typename Terms, typename Inner>
typename Scale::PlainObject canonical_row_floor(const Matrix& j, const Scale& physical_scale,
    const Terms& terms, const Inner& inner) {
    const double epsilon = std::numeric_limits<double>::epsilon();
    const double gamma = j.cols() * epsilon / (1 - j.cols() * epsilon);
    return (gamma * (j.cwiseAbs() * physical_scale)).cwiseMax(inner + epsilon * terms);
}

template <typename Vector>
struct LinearError final {
    Vector correction, floor, defect, radius;
    std::size_t rank = 0;
    double threshold = 0, matrix_norm = 0, inverse_norm = INFINITY, condition = INFINITY, backward = INFINITY;
    std::string cause;
};

template <typename Jacobian, typename Residual, typename Projection, typename Floor, typename RowScale, typename ColumnScale>
LinearError<typename Residual::PlainObject> linear_error(const Jacobian& j, const Residual& f, const Projection& g,
    const Floor& eta, const RowScale& row, const ColumnScale& column) {
    using Matrix = typename Jacobian::PlainObject;
    using Vector = typename Residual::PlainObject;
    LinearError<Vector> out;
    const Matrix a = row.asDiagonal() * j * column.asDiagonal();
    if (a.rows() != a.cols() || !a.allFinite() || !f.allFinite() || !eta.allFinite()
        || !g.allFinite() || (eta.array() < 0).any()) { out.cause = "nonfinite_or_unavailable_rows"; return out; }
    const auto lu = a.fullPivLu();
    out.rank = static_cast<std::size_t>(lu.rank()); out.threshold = lu.threshold();
    out.matrix_norm = a.cwiseAbs().rowwise().sum().maxCoeff();
    if (out.rank != static_cast<std::size_t>(a.cols())) { out.cause = "rank_deficient"; return out; }
    const Vector right = -row.cwiseProduct(f), d = lu.solve(right);
    const Vector defect = compensated_defect(a, d, right);
    const Matrix gc = g * column.asDiagonal();
    out.correction = gc * d;
    out.floor = Vector::Zero(g.rows()); out.defect = out.floor;
    Vector sums = Vector::Zero(a.cols());
    for (std::ptrdiff_t i = 0; i < a.cols(); ++i) {
        const Vector unit = lu.solve(Vector::Unit(a.cols(), i));
        sums += unit.cwiseAbs();
        const Vector action = (gc * unit).cwiseAbs();
        out.floor += action * (std::abs(row[i]) * eta[i]);
        out.defect += action * std::abs(defect[i]);
    }
    out.inverse_norm = sums.maxCoeff(); out.condition = out.matrix_norm * out.inverse_norm;
    out.backward = backward_error(defect.cwiseAbs().maxCoeff(), out.matrix_norm, d.cwiseAbs().maxCoeff(), right.cwiseAbs().maxCoeff());
    out.radius = out.correction.cwiseAbs() + out.floor + out.defect;
    if (!out.radius.allFinite()) out.cause = "nonfinite_projection";
    return out;
}

enum class Quantity { Linear, Positive, LogPositive };
enum class ComparisonMode { Strict, EstimatedError };
struct NumericalEstimate final {
    double absolute = NAN, chart = NAN;
    bool available = false, resolved = false;
    std::string cause = "not_evaluated";
};
struct ComparisonEligibility final {
    bool declaration = true, branch = true, phase_matching = true, inventory = true;
    bool residual_a = true, residual_b = true, admission_a = true, admission_b = true;
    bool valid() const { return declaration && branch && phase_matching && inventory && residual_a && residual_b && admission_a && admission_b; }
};
struct EquilibriumComparison final {
    std::string judgment = "failed", cause;
    double difference = 0, chart_difference = 0, radius_sum = 0;
};
inline double log_difference(double a, double b) {
    const double relative = (a - b) / b;
    return std::abs(std::isfinite(relative) && std::abs(relative) < policy::LOG_RATIO_NEAR_UNITY.value ? std::log1p(relative) : std::log(a) - std::log(b));
}
inline NumericalEstimate project_numerical_estimate(double value, double mapped, double radius, bool positive, bool constant) {
    core::Scope scope("project_numerical_estimate");
    const double rounding = constant ? 0.0 : std::numeric_limits<double>::epsilon() * (positive ? 1.0 : std::abs(value));
    const double displacement = positive ? log_difference(value, mapped) : std::abs(value - mapped);
    const double chart = radius + displacement + rounding;
    const double absolute = positive ? value * std::expm1(chart) : chart;
    const bool available = std::isfinite(absolute) && std::isfinite(chart);
    return {absolute, chart, available, available, available ? "" : "nonfinite_projection"};
}
inline EquilibriumComparison compare_equilibrium_quantity(Quantity quantity, double a, double b, double allowance,
    const NumericalEstimate& ea, const NumericalEstimate& eb, ComparisonEligibility eligibility, ComparisonMode mode) {
    core::Scope scope("compare_equilibrium_quantity");
    EquilibriumComparison out; out.difference = std::abs(a - b);
    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(allowance) || allowance < 0) { out.cause = "nonfinite_operand"; return out; }
    if (quantity == Quantity::Positive && (!(a > 0) || !(b > 0))) { out.cause = "positive_domain"; return out; }
    out.chart_difference = quantity == Quantity::Positive ? log_difference(a, b) : out.difference;
    if (!eligibility.valid()) { out.cause = "prerequisite"; return out; }
    if (out.difference <= allowance) { out.judgment = "ordinary"; out.cause = "ordinary_allowance"; return out; }
    if (mode == ComparisonMode::Strict) { out.cause = "strict_comparison"; return out; }
    const auto usable = [](const NumericalEstimate& e) {
        return e.available && e.resolved && std::isfinite(e.absolute) && std::isfinite(e.chart) && e.absolute >= 0 && e.chart >= 0;
    };
    if (!usable(ea) || !usable(eb)) { out.cause = "unavailable_or_unresolved"; return out; }
    if (ea.absolute <= allowance && eb.absolute <= allowance) { out.cause = "ordinary_achievable"; return out; }
    out.radius_sum = quantity == Quantity::Linear ? ea.absolute + eb.absolute : ea.chart + eb.chart;
    if (!std::isfinite(out.radius_sum)) { out.cause = "nonfinite_sum"; return out; }
    if (out.chart_difference <= out.radius_sum) {
        out.judgment = "estimated_error"; out.cause = "consistent_within_estimated_numerical_error";
    } else out.cause = "outside_estimates";
    return out;
}
}  // namespace epcsaft::core
