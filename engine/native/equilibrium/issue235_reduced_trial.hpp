#pragma once

#include "equilibrium/gibbs_minimizer.hpp"
#include "core/issue235_prototype.hpp"

#include <Eigen/Dense>

namespace epcsaft::equilibrium::issue235_trial {

class Chart final : public Ipopt::TNLP {
    using Index = Ipopt::Index;
    using Number = Ipopt::Number;
    using Matrix = Eigen::MatrixXd;
    using Vector = Eigen::VectorXd;
    using RowMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

public:
    Chart(GibbsMinimizer& full, const CompiledSystem& c) : full_(full) {
        Index jac, h;
        IndexStyleEnum style;
        full.get_nlp_info(full_n_, full_m_, jac, h, style);
        const auto& support = c.problem.phases[0].support;
        composition_size_ = support.size();
        if (full_n_ != composition_size_ + 2 || !c.shape[0].amount ||
            *c.shape[0].amount != static_cast<std::size_t>(composition_size_) ||
            c.shape[0].volume != static_cast<std::size_t>(composition_size_ + 1)) return;
        for (Index i = 0; i < composition_size_; ++i)
            if (c.shape[0].x[i] != static_cast<std::size_t>(i)) return;

        const auto& charge = c.model->parameters().charge_number;
        const bool charged = std::any_of(support.begin(), support.end(), [&](auto i) { return charge[i] != 0; });
        const Index rank = charged ? 2 : 1;
        if (composition_size_ <= rank || full_m_ != rank + 1) return;
        if (std::count_if(c.rows.begin(), c.rows.end(), [](const auto& r) { return r.kind == RowKind::Normalization; }) != 1) return;
        for (const auto& r : c.rows)
            if (!stationarity_row(r.kind) && r.kind != RowKind::Normalization
                && r.kind != RowKind::Moiety && r.kind != RowKind::Electro) return;

        g_.resize(full_m_);
        full_lower_.resize(full_n_); full_upper_.resize(full_n_);
        point_.resize(full_n_); weights_.resize(full_m_);
        if (!full.get_bounds_info(full_n_, full_lower_.data(), full_upper_.data(), full_m_, g_.data(), weights_.data())
            || !full.get_starting_point(full_n_, true, point_.data(), false, nullptr, nullptr, full_m_, false, nullptr))
            return;
        // Optimized stability starts have unit amount in both the original and TNLP charts.
        if (point_[composition_size_] != 1.0 || std::abs(c.feed_total - 1.0) > roundoff) return;
        amount_ = c.feed_total;
        volume_start_ = point_[composition_size_ + 1];

        Matrix inventory(rank, composition_size_);
        inventory.row(0).setOnes();
        if (charged) for (Index i = 0; i < composition_size_; ++i) inventory(1, i) = charge[support[i]];
        const Vector seed = point_.head(composition_size_).array().exp();
        const auto factor = (inventory * seed.asDiagonal()).colPivHouseholderQr();
        if (factor.rank() != rank) return;
        const auto permutation = factor.colsPermutation().indices();
        Matrix pivots(rank, rank);
        for (Index i = 0; i < rank; ++i) pivots.col(i) = inventory.col(permutation[i]);
        const auto pivot = pivots.fullPivLu();
        if (pivot.rank() != rank || pivot.rcond() <= roundoff) return;
        Vector target = Vector::Zero(rank);
        target[0] = 1.0;
        const Vector correction = pivot.solve(target - inventory * seed);
        origin_ = seed;
        for (Index i = 0; i < rank; ++i) origin_[permutation[i]] += correction[i];
        if ((origin_ - seed).cwiseAbs().maxCoeff() > roundoff || (origin_.array() <= 0).any()) return;

        const Index free = composition_size_ - rank;
        basis_ = Matrix::Zero(composition_size_, free);
        for (Index j = 0; j < free; ++j) {
            const auto i = permutation[rank + j];
            const Vector column = pivot.solve(-inventory.col(i));
            basis_(i, j) = 1.0;
            for (Index k = 0; k < rank; ++k) basis_(permutation[k], j) = column[k];
        }
        if ((inventory * basis_).cwiseAbs().maxCoeff() > roundoff
            || (inventory * origin_ - target).cwiseAbs().maxCoeff() > roundoff) return;
        // Check all actual balance equations, including the compiled feed admission.
        Matrix equations(full_m_ - 1, composition_size_);
        Index row = 0;
        for (const auto& r : c.rows) {
            if (r.kind != RowKind::Moiety && r.kind != RowKind::Electro) continue;
            for (Index i = 0; i < composition_size_; ++i)
                equations(row, i) = r.kind == RowKind::Moiety ? c.B[r.index][support[i]] : charge[support[i]];
            const double rhs = r.kind == RowKind::Moiety ? c.balance[r.index] : 0.0;
            const double scale = std::max(1.0, equations.row(row).cwiseAbs().maxCoeff());
            if (std::abs(equations.row(row).dot(origin_) - rhs / amount_) > roundoff * scale
                || (equations.row(row) * basis_).cwiseAbs().maxCoeff() > roundoff * scale) return;
            ++row;
        }
        if (row != full_m_ - 1) return;
        Matrix independence(full_m_, composition_size_ + 1);
        independence.topLeftCorner(row, composition_size_) = equations;
        independence.topRightCorner(row, 1).setZero();
        row = 0;
        for (const auto& r : c.rows) {
            if (r.kind == RowKind::Moiety) independence(row++, composition_size_) = c.balance[r.index] / amount_;
            if (r.kind == RowKind::Electro) ++row;
        }
        independence.bottomRows(1).setOnes();
        independence(full_m_ - 1, composition_size_) = 0.0;
        if (independence.fullPivLu().rank() != full_m_) return;

        n_ = free + 1;
        jacobian_ = Matrix::Zero(full_n_, n_);
        jacobian_(composition_size_ + 1, n_ - 1) = 1.0;
        gradient_.resize(full_n_);
        hessian_.resize(full_n_, full_n_);
        packed_.resize(full_n_ * (full_n_ + 1) / 2);
        full_jacobian_.resize(full_m_, full_n_);
        admitted = (origin_.array() >= full_lower_.head(composition_size_).array().exp()).all()
            && (origin_.array() <= full_upper_.head(composition_size_).array().exp()).all()
            && volume_start_ >= full_lower_[full_n_ - 1] && volume_start_ <= full_upper_[full_n_ - 1];
    }

    bool get_nlp_info(Index& n, Index& m, Index& jac, Index& h, IndexStyleEnum& style) override {
        n = n_; m = composition_size_; jac = m * (n - 1); h = n * (n + 1) / 2;
        style = C_STYLE;
        return admitted;
    }

    bool get_bounds_info(Index n, Number* lower, Number* upper, Index m, Number* gl, Number* gu) override {
        std::fill_n(lower, n, -1.0e19);
        std::fill_n(upper, n, 1.0e19);
        lower[n - 1] = full_lower_[full_n_ - 1]; upper[n - 1] = full_upper_[full_n_ - 1];
        for (Index i = 0; i < m; ++i) {
            gl[i] = std::exp(full_lower_[i]);
            gu[i] = std::exp(full_upper_[i]);
        }
        return true;
    }

    bool get_starting_point(Index n, bool ix, Number* x, bool iz, Number*, Number*,
                            Index, bool il, Number*) override {
        if (!ix || iz || il) return false;
        std::fill_n(x, n - 1, 0.0);
        x[n - 1] = volume_start_;
        return true;
    }

    bool eval_f(Index, const Number* x, bool, Number& value) override {
        return lift(x) && full_.eval_f(full_n_, point_.data(), true, value);
    }

    bool eval_grad_f(Index, const Number* x, bool, Number* out) override {
        if (!lift(x) || !gradient()) return false;
        Eigen::Map<Vector>(out, n_) = jacobian_.transpose() * gradient_;
        return true;
    }

    bool eval_g(Index, const Number* x, bool, Index m, Number* g) override {
        Eigen::Map<Vector>(g, m) = composition(x);
        return true;
    }

    bool eval_jac_g(Index n, const Number*, bool, Index m, Index, Index* rows,
                    Index* cols, Number* values) override {
        Index k = 0;
        for (Index i = 0; i < m; ++i) for (Index j = 0; j < n - 1; ++j) {
            if (values) values[k] = basis_(i, j);
            else { rows[k] = i; cols[k] = j; }
            ++k;
        }
        return true;
    }

    bool eval_h(Index n, const Number* x, bool, Number factor, Index, const Number*,
                bool, Index, Index* rows, Index* cols, Number* values) override {
        Matrix reduced;
        if (values) {
            weights_.setZero();
            if (!lift(x) || !gradient() || !full_.eval_h(full_n_, point_.data(), true, factor,
                    full_m_, weights_.data(), true, packed_.size(), nullptr, nullptr, packed_.data())) return false;
            Index k = 0;
            for (Index i = 0; i < full_n_; ++i) for (Index j = 0; j <= i; ++j)
                hessian_(i, j) = packed_[k++];
            reduced = jacobian_.transpose() * hessian_.selfadjointView<Eigen::Lower>() * jacobian_;
            for (Index i = 0; i < composition_size_; ++i)
                reduced.noalias() -= (factor * gradient_[i]) * jacobian_.row(i).transpose() * jacobian_.row(i);
        }
        Index k = 0;
        for (Index i = 0; i < n; ++i) for (Index j = 0; j <= i; ++j) {
            if (values) values[k] = reduced(i, j);
            else { rows[k] = i; cols[k] = j; }
            ++k;
        }
        return true;
    }

    bool intermediate_callback(Ipopt::AlgorithmMode, Index iteration, Number, Number, Number, Number,
                               Number, Number, Number, Number, Index, const Ipopt::IpoptData*,
                               Ipopt::IpoptCalculatedQuantities*) override {
        iterations = iteration;
        return true;
    }

    void finalize_solution(Ipopt::SolverReturn status, Index n, const Number* x, const Number* zl,
        const Number* zu, Index m, const Number*, const Number* lambda, Number objective,
        const Ipopt::IpoptData* data, Ipopt::IpoptCalculatedQuantities* calculated) override {
        if (status != Ipopt::SUCCESS || !lift(x) || !gradient()) return;
        if ((point_.array() < full_lower_.array()).any() || (point_.array() > full_upper_.array()).any()) return;
        if (!full_.eval_jac_g(full_n_, point_.data(), true, full_m_, full_n_ * full_m_,
                nullptr, nullptr, full_jacobian_.data())) return;
        Vector lower = Vector::Zero(full_n_), upper = lower;
        lower[full_n_ - 1] = zl[n - 1]; upper[full_n_ - 1] = zu[n - 1];
        for (Index i = 0; i < m; ++i) {
            lower[i] = std::max(0.0, -lambda[i]) * fractions_[i];
            upper[i] = std::max(0.0, lambda[i]) * fractions_[i];
        }
        const auto factor = full_jacobian_.transpose().colPivHouseholderQr();
        if (factor.rank() != full_m_) return;
        // The native full-driver postprocessor applies its usual multiplier lift and endpoint checks.
        const Vector equality = factor.solve(-gradient_ + lower - upper);
        full_.finalize_solution(status, full_n_, point_.data(), lower.data(), upper.data(),
                               full_m_, nullptr, equality.data(), objective, data, calculated);
    }

    bool admitted = false;
    std::size_t iterations = 0;

private:
    static constexpr double roundoff = 64 * std::numeric_limits<double>::epsilon();
    Vector composition(const Number* x) const {
        return origin_ + basis_ * Eigen::Map<const Vector>(x, n_ - 1);
    }

    bool lift(const Number* x) {
        fractions_ = composition(x);
        if ((fractions_.array() <= 0.0).any() || !fractions_.allFinite()) return false;
        point_.head(composition_size_) = fractions_.array().log();
        point_[composition_size_] = amount_;
        point_[full_n_ - 1] = x[n_ - 1];
        jacobian_.topLeftCorner(composition_size_, n_ - 1) = fractions_.cwiseInverse().asDiagonal() * basis_;
        return true;
    }

    bool gradient() { return full_.eval_grad_f(full_n_, point_.data(), true, gradient_.data()); }
    GibbsMinimizer& full_;
    Index full_n_ = 0, full_m_ = 0, n_ = 0, composition_size_ = 0;
    double amount_ = 0.0, volume_start_ = 0.0;
    Vector point_, fractions_, origin_, full_lower_, full_upper_, gradient_, packed_, weights_, g_;
    Matrix basis_, jacobian_, hessian_;
    RowMatrix full_jacobian_;
};

}
