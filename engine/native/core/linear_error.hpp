#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <string>

namespace epcsaft::core {

// Rigal-Gaches normwise backward error; the zero right-hand side needs no unit floor.
inline double backward_error(double residual, double jacobian_norm, double solution, double right) {
    const double denominator = jacobian_norm * solution + right;
    return residual == 0.0 ? 0.0 : denominator > 0.0 ? residual / denominator : std::numeric_limits<double>::infinity();
}

// Ogita-Rump-Oishi Dot2: FMA TwoProduct and TwoSum preserve the cancellation defect.
inline Eigen::VectorXd compensated_defect(const Eigen::MatrixXd& matrix, const Eigen::VectorXd& solution, const Eigen::VectorXd& right) {
    Eigen::VectorXd result(right.size());
    for (Eigen::Index row = 0; row < right.size(); ++row) {
        double sum = right(row), error = 0.0;
        for (Eigen::Index column = 0; column < solution.size(); ++column) {
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
    for (Eigen::Index i = 0; i < columns.size(); ++i) {
        const Vector rhs = Vector::Unit(columns.size(), i);
        sums += (columns.asDiagonal() * lu.solve(rows.asDiagonal() * rhs)).cwiseAbs();
    }
    return sums.maxCoeff();
}

// Propagate nonnegative row errors through an already admitted factorization.
template <typename Factor>
Eigen::VectorXd inverse_error_action(const Factor& lu, const Eigen::VectorXd& right) {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(right.size());
    for (Eigen::Index i = 0; i < right.size(); ++i)
        result += lu.solve(Eigen::VectorXd::Unit(right.size(), i)).cwiseAbs() * right[i];
    return result;
}

template <typename Real>
Real reduced_error_envelope(Real inverse, Real residual, Real floor, Real displacement = 0) {
    return inverse * (residual + floor) + displacement;
}

inline Eigen::VectorXd canonical_row_floor(const Eigen::MatrixXd& j, const Eigen::VectorXd& physical_scale,
    const Eigen::VectorXd& terms, const Eigen::VectorXd& inner) {
    const double epsilon = std::numeric_limits<double>::epsilon();
    const double gamma = j.cols() * epsilon / (1 - j.cols() * epsilon);
    return (gamma * (j.cwiseAbs() * physical_scale)).cwiseMax(inner + epsilon * terms);
}

struct LinearError final {
    Eigen::VectorXd correction, floor, defect, radius;
    std::size_t rank = 0;
    double threshold = 0, matrix_norm = 0, inverse_norm = INFINITY, condition = INFINITY, backward = INFINITY;
    std::string cause;
};

inline LinearError linear_error(const Eigen::MatrixXd& j, const Eigen::VectorXd& f, const Eigen::MatrixXd& g,
    const Eigen::VectorXd& eta, const Eigen::VectorXd& row, const Eigen::VectorXd& column) {
    LinearError out;
    const Eigen::MatrixXd a = row.asDiagonal() * j * column.asDiagonal();
    if (a.rows() != a.cols() || !a.allFinite() || !f.allFinite() || !eta.allFinite()
        || !g.allFinite() || (eta.array() < 0).any()) { out.cause = "nonfinite_or_unavailable_rows"; return out; }
    const Eigen::FullPivLU<Eigen::MatrixXd> lu(a);
    out.rank = static_cast<std::size_t>(lu.rank()); out.threshold = lu.threshold();
    out.matrix_norm = a.cwiseAbs().rowwise().sum().maxCoeff();
    if (out.rank != static_cast<std::size_t>(a.cols())) { out.cause = "rank_deficient"; return out; }
    const Eigen::VectorXd right = -row.cwiseProduct(f), d = lu.solve(right);
    const Eigen::VectorXd defect = compensated_defect(a, d, right);
    const Eigen::MatrixXd gc = g * column.asDiagonal();
    out.correction = gc * d;
    out.floor = Eigen::VectorXd::Zero(g.rows()); out.defect = out.floor;
    Eigen::VectorXd sums = Eigen::VectorXd::Zero(a.cols());
    for (Eigen::Index i = 0; i < a.cols(); ++i) {
        const Eigen::VectorXd unit = lu.solve(Eigen::VectorXd::Unit(a.cols(), i));
        sums += unit.cwiseAbs();
        const Eigen::VectorXd action = (gc * unit).cwiseAbs();
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
    double absolute = 0, chart = 0;
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
    return std::abs(std::isfinite(relative) && std::abs(relative) < 0.5 ? std::log1p(relative) : std::log(a) - std::log(b));
}
inline EquilibriumComparison compare_equilibrium_quantity(Quantity quantity, double a, double b, double allowance,
    const NumericalEstimate& ea, const NumericalEstimate& eb, ComparisonEligibility eligibility, ComparisonMode mode) {
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
