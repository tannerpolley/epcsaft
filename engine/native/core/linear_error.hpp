#pragma once
#include "core/wider_scalar.hpp"

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

// Eigen's supplied-factor estimator observes the 1-norm after rank admission.
template <typename Factor>
NumericalRank conditioning(const Factor& factor, const char* basis, std::size_t iteration = 0,
    const char* point_kind = "last_linearization", const std::vector<double>& controls = {}) {
    auto record = numerical_rank(factor, basis);
    ConditioningObservation observation{basis, "", "1", "Eigen supplied-factor estimate", "not_evaluated", point_kind};
    observation.iteration = iteration; observation.controls = controls;
    observation.precision_bits = std::numeric_limits<typename Factor::MatrixType::Scalar>::digits;
    if (factor.rows() == 0 || factor.rows() != factor.cols() || !factor.matrixLU().allFinite())
        observation.cause = "nonfinite_or_unavailable_rows";
    else if (record.rank != record.dimension) observation.cause = "rank_deficient";
    else {
        const auto reciprocal = factor.rcond();
        observation.condition = double(1 / reciprocal);
        observation.available = std::isfinite(observation.condition) && reciprocal > 0;
        observation.cause = observation.available ? "physical_error_operands_unavailable" : "nonfinite_inverse";
    }
    observe_conditioning(record, std::move(observation));
    return record;
}

// Rigal-Gaches normwise backward error; the zero right-hand side needs no unit floor.
template <typename Real>
Real backward_error(Real residual, Real jacobian_norm, Real solution, Real right) {
    const Real denominator = jacobian_norm * solution + right;
    return residual == 0.0 ? 0.0 : denominator > 0.0 ? residual / denominator : std::numeric_limits<double>::infinity();
}

#ifdef EPCSAFT_SOFTWARE_WIDER
inline WiderScalar product_residual(WiderScalar a, WiderScalar b, WiderScalar p) {
    using Product = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<128,
        boost::multiprecision::digit_base_2, void, int, -32766, 32767>, boost::multiprecision::et_off>;
    const Product residual = Product(a) * Product(b) - Product(p);
    if (!wide_math::isfinite(residual) || wide_math::abs(residual) > Product(std::numeric_limits<WiderScalar>::max()))
        throw std::runtime_error("product residual outside wider scalar range");
    return static_cast<WiderScalar>(residual);
}
#endif

// Ogita-Rump-Oishi Dot2: FMA TwoProduct and TwoSum preserve the cancellation defect.
template <typename Matrix, typename Solution, typename Right>
typename Right::PlainObject compensated_defect(const Matrix& matrix, const Solution& solution, const Right& right) {
    typename Right::PlainObject result(right.size());
    for (std::ptrdiff_t row = 0; row < right.size(); ++row) {
        using Real = typename Right::Scalar;
        Real sum = right(row), error = 0.0;
        for (std::ptrdiff_t column = 0; column < solution.size(); ++column) {
            const Real term = -matrix(row, column) * solution(column), next = sum + term, part = next - sum;
#ifdef EPCSAFT_SOFTWARE_WIDER
            if constexpr (std::is_same_v<Real, WiderScalar>)
                error += (sum - (next - part)) + (term - part) + product_residual(-matrix(row, column), solution(column), term);
            else
#endif
            error += (sum - (next - part)) + (term - part) + std::fma(-matrix(row, column), solution(column), -term);
            sum = next;
        }
        result(row) = sum + error;
    }
    return result;
}

template <typename Matrix, typename Scale, typename Terms, typename Inner>
typename Scale::PlainObject canonical_row_floor(const Matrix& j, const Scale& physical_scale,
    const Terms& terms, const Inner& inner) {
    using Real = typename Matrix::Scalar;
    const Real epsilon = std::numeric_limits<Real>::epsilon();
    const Real gamma = j.cols() * epsilon / (1 - j.cols() * epsilon);
    return (gamma * (j.cwiseAbs() * physical_scale)).cwiseMax(inner + epsilon * terms);
}

template <typename Vector>
struct LinearError final {
    Vector correction, floor, defect, radius, row_floor;
    std::size_t rank = 0;
    typename Vector::Scalar threshold = 0, matrix_norm = 0, inverse_norm = INFINITY, condition = INFINITY, backward = INFINITY;
    std::string cause;
};

template <typename Jacobian, typename Residual, typename Projection, typename Floor, typename RowScale, typename ColumnScale,
    typename Factor = std::decay_t<decltype(std::declval<typename Jacobian::PlainObject>().fullPivLu())>>
LinearError<typename Residual::PlainObject> linear_error(const Jacobian& j, const Residual& f, const Projection& g,
    const Floor& eta, const RowScale& row, const ColumnScale& column,
    const Factor* supplied = nullptr,
    const typename Residual::PlainObject* correction = nullptr) {
    using Matrix = typename Jacobian::PlainObject;
    using Vector = typename Residual::PlainObject;
    LinearError<Vector> out; out.row_floor = eta;
    if (j.rows() == 0 || j.rows() != j.cols() || f.size() != j.rows() || eta.size() != j.rows()
        || row.size() != j.rows() || column.size() != j.cols() || g.cols() != j.cols()
        || !row.allFinite() || !column.allFinite() || (correction && (correction->size() != j.cols() || !correction->allFinite()))) {
        out.cause = "nonfinite_or_unavailable_rows"; return out;
    }
    const Matrix a = row.asDiagonal() * j * column.asDiagonal();
    if (a.rows() != a.cols() || !a.allFinite() || !f.allFinite() || !eta.allFinite()
        || !g.allFinite() || (eta.array() < 0).any()) { out.cause = "nonfinite_or_unavailable_rows"; return out; }
    std::optional<Factor> owned;
    const auto& lu = supplied ? *supplied : owned.emplace(a);
    if (lu.rows() != a.rows() || lu.cols() != a.cols()) { out.cause = "factor_shape"; return out; }
    out.rank = static_cast<std::size_t>(lu.rank()); out.threshold = lu.threshold();
    out.matrix_norm = a.cwiseAbs().rowwise().sum().maxCoeff();
    if (out.rank != static_cast<std::size_t>(a.cols())) { out.cause = "rank_deficient"; return out; }
    const Vector right = -row.cwiseProduct(f), d = correction ? *correction : Vector(lu.solve(right));
    const Vector defect = compensated_defect(a, d, right);
    const Matrix gc = g * column.asDiagonal();
    out.correction = gc * d;
    out.floor = Vector::Zero(g.rows()); out.defect = out.floor;
    Vector sums = Vector::Zero(a.cols());
    for (std::ptrdiff_t i = 0; i < a.cols(); ++i) {
        const Vector unit = lu.solve(Vector::Unit(a.cols(), i));
        sums += unit.cwiseAbs();
        const Vector action = (gc * unit).cwiseAbs();
        out.floor += action * (wide_math::abs(row[i]) * eta[i]);
        out.defect += action * wide_math::abs(defect[i]);
    }
    out.inverse_norm = sums.maxCoeff(); out.condition = out.matrix_norm * out.inverse_norm;
    out.backward = backward_error(defect.cwiseAbs().maxCoeff(), out.matrix_norm, d.cwiseAbs().maxCoeff(), right.cwiseAbs().maxCoeff());
    out.radius = out.correction.cwiseAbs() + out.floor + out.defect;
    if (!out.radius.allFinite()) out.cause = "nonfinite_projection";
    return out;
}

// Preserve the polish admission operand: infinity-norm estimation uses the transpose.
template <typename Matrix, typename Residual, typename Correction, typename Scale>
double polish_forward_bound(const Matrix& j, const Residual& r, const Correction& d, const Scale& scale) {
    const double unavailable = std::numeric_limits<double>::infinity();
    if (j.rows() == 0 || j.rows() != j.cols() || r.size() != j.rows() || d.size() != j.cols()
        || scale.size() != j.cols() || !j.allFinite() || !d.allFinite() || !r.allFinite()
        || !scale.allFinite() || !(scale.minCoeff() > 0.0)) return unavailable;
    const typename Matrix::PlainObject scaled = j * scale.asDiagonal();
    const auto lu = scaled.transpose().eval().fullPivLu();
    if (lu.rank() != j.cols()) return unavailable;
    const double denominator = lu.rcond() * scaled.cwiseAbs().rowwise().sum().maxCoeff();
    if (!(denominator > 0.0) || !std::isfinite(denominator)) return unavailable;
    return d.cwiseQuotient(scale).cwiseAbs().maxCoeff()
        + (j * d + r).cwiseAbs().maxCoeff() / denominator;
}

}  // namespace epcsaft::core
