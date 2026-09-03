#pragma once

#include "equilibrium/problem.hpp"

#include <IpTNLP.hpp>

#include <cstddef>
#include <optional>
#include <vector>

namespace epcsaft::core { class Scope; }

namespace epcsaft::equilibrium {

struct NlpShape final {
    Ipopt::Index variables = 0;
    Ipopt::Index constraints = 0;
    Ipopt::Index jacobian_entries = 0;
    Ipopt::Index hessian_entries = 0;
};

std::vector<double> boundary_start(const CompiledSystem& compiled);

std::vector<double> solve_implicit(
    const CompiledSystem& compiled,
    const std::vector<double>& point,
    const std::vector<double>& parameter_rows
);

std::vector<double> sensitivity(
    const CompiledSystem& compiled,
    const std::vector<double>& point,
    const std::vector<double>& lever
);

struct SensitivityInfo final {
    std::size_t equations = 0;
    std::size_t parameters = 0;
    std::size_t rank = 0;
    double condition_estimate = 0.0;
    bool nonsingular = false;
};

SensitivityInfo inspect_sensitivity(
    const CompiledSystem& compiled,
    const std::vector<double>& point,
    const std::vector<double>& parameter_rows
);

class GibbsMinimizer final : public Ipopt::TNLP {
public:
    GibbsMinimizer(const CompiledSystem& compiled, std::vector<double> start,
                   core::Scope* scope = nullptr);
    bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m,
                      Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
                      IndexStyleEnum& index_style) override;
    bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
                         Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u) override;
    bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
                            bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
                            Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda) override;
    bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                Ipopt::Number& value) override;
    bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                     Ipopt::Number* grad_f) override;
    bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                Ipopt::Index m, Ipopt::Number* g) override;
    bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                    Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
                    Ipopt::Index* jCol, Ipopt::Number* values) override;
    bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
                bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                Ipopt::Index* jCol, Ipopt::Number* values) override;
    bool intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter,
                                Ipopt::Number obj_value, Ipopt::Number inf_pr,
                                Ipopt::Number inf_du, Ipopt::Number mu,
                                Ipopt::Number d_norm, Ipopt::Number regularization_size,
                                Ipopt::Number alpha_du, Ipopt::Number alpha_pr,
                                Ipopt::Index ls_trials, const Ipopt::IpoptData* ip_data,
                                Ipopt::IpoptCalculatedQuantities* ip_cq) override;
    void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n,
                           const Ipopt::Number* x, const Ipopt::Number* z_L,
                           const Ipopt::Number* z_U, Ipopt::Index m,
                           const Ipopt::Number* g, const Ipopt::Number* lambda,
                           Ipopt::Number objective,
                           const Ipopt::IpoptData* ip_data,
                           Ipopt::IpoptCalculatedQuantities* ip_cq) override;

    const std::vector<double>& solution() const noexcept { return solution_; }
    const std::vector<double>& constraint_multipliers() const noexcept {
        return constraint_multipliers_;
    }
    Ipopt::SolverReturn status() const noexcept { return status_; }
    std::size_t iterations() const noexcept { return iterations_; }
    std::size_t evaluations() const noexcept { return evaluations_; }
    std::size_t hessian_calls() const noexcept { return hessian_calls_; }

    double regularization_size() const noexcept { return regularization_size_; }

    NlpShape shape() const {
        const Ipopt::Index variables = static_cast<Ipopt::Index>(
            compiled_.pinned_count == 0
                ? compiled_.primal_size
                : compiled_.primal_size + compiled_.rank + compiled_.charge_row_count);
        const Ipopt::Index constraints = static_cast<Ipopt::Index>(constraint_rows().size());
        return {variables, constraints, variables * constraints,
                variables * (variables + 1) / 2};
    }

    bool uses_exact_hessian() const noexcept { return true; }

    bool has_zero_objective() const noexcept { return compiled_.pinned_count != 0; }

    const CompiledSystem& compiled_system() const noexcept { return compiled_; }

    std::size_t primal_size() const noexcept { return compiled_.primal_size; }

    std::size_t active_row_count() const noexcept { return constraint_rows().size(); }

    std::size_t multiplier_count() const noexcept {
        return compiled_.rank + compiled_.charge_row_count;
    }

    bool is_square() const noexcept { return compiled_.is_square(); }

    bool topology_is_declared() const noexcept {
        return compiled_.problem.topology_declared;
    }

    const std::vector<double>& starting_point() const noexcept { return start_; }

private:
    std::vector<double> point(Ipopt::Index n, const Ipopt::Number* x) const;
    std::vector<std::size_t> constraint_rows() const;
    std::vector<double> objective_gradient(const std::vector<double>& x) const;

    const CompiledSystem& compiled_;
    std::vector<double> start_;
    std::vector<double> solution_;
    std::vector<double> constraint_multipliers_;
    Ipopt::SolverReturn status_ = Ipopt::UNASSIGNED;
    std::size_t iterations_ = 0;
    std::size_t evaluations_ = 0;
    std::size_t hessian_calls_ = 0;
    core::Scope* scope_ = nullptr;
    double regularization_size_ = 0.0;
};

}  // namespace epcsaft::equilibrium
