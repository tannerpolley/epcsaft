#include "chemical_equilibrium.hpp"
#include "provider.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpTNLP.hpp>

namespace epcsaft_engine::equilibrium {
namespace {

constexpr double kGasConstantJPerMolK = 8.31446261815324;

std::string ipopt_status_name(Ipopt::ApplicationReturnStatus status) {
    switch (status) {
        case Ipopt::Solve_Succeeded:
            return "solve_succeeded";
        case Ipopt::Solved_To_Acceptable_Level:
            return "solved_to_acceptable_level";
        case Ipopt::Infeasible_Problem_Detected:
            return "infeasible_problem_detected";
        case Ipopt::Search_Direction_Becomes_Too_Small:
            return "search_direction_too_small";
        case Ipopt::Diverging_Iterates:
            return "diverging_iterates";
        case Ipopt::Maximum_Iterations_Exceeded:
            return "maximum_iterations_exceeded";
        case Ipopt::Restoration_Failed:
            return "restoration_failed";
        case Ipopt::Error_In_Step_Computation:
            return "error_in_step_computation";
        case Ipopt::Invalid_Problem_Definition:
            return "invalid_problem_definition";
        case Ipopt::Invalid_Number_Detected:
            return "invalid_number_detected";
        case Ipopt::Unrecoverable_Exception:
            return "unrecoverable_exception";
        default:
            return "ipopt_status_" + std::to_string(static_cast<int>(status));
    }
}

double inf_norm(const std::vector<double>& values) {
    double result = 0.0;
    for (double value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

struct NormalizedBalances {
    DenseMatrix matrix;
    std::vector<double> totals;
    std::vector<double> scales;
};

NormalizedBalances normalized_balances(const CompiledReactionSystem& system) {
    NormalizedBalances result{
        system.balance_matrix,
        system.balance_totals,
        std::vector<double>(system.balance_matrix.rows, 1.0),
    };
    for (std::size_t row = 0; row < result.matrix.rows; ++row) {
        double scale = std::abs(result.totals[row]);
        for (std::size_t species = 0; species < result.matrix.columns; ++species) {
            scale += std::abs(
                result.matrix(row, species) * system.feed_amounts[species]
            );
        }
        scale = std::max(scale, std::numeric_limits<double>::min());
        result.scales[row] = scale;
        result.totals[row] /= scale;
        for (std::size_t species = 0; species < result.matrix.columns; ++species) {
            result.matrix(row, species) /= scale;
        }
    }
    return result;
}

std::size_t declared_direction_rank(const CompiledReactionSystem& system) {
    const std::size_t species_count = system.species_count;
    const std::size_t reaction_count = system.reaction_matrix.rows;
    Eigen::MatrixXd directions = Eigen::MatrixXd::Zero(
        static_cast<Eigen::Index>(2 * reaction_count + species_count),
        static_cast<Eigen::Index>(2 * species_count)
    );
    for (std::size_t reaction = 0; reaction < reaction_count; ++reaction) {
        for (std::size_t species = 0; species < species_count; ++species) {
            const double coefficient = system.reaction_matrix(reaction, species);
            directions(
                static_cast<Eigen::Index>(reaction),
                static_cast<Eigen::Index>(species)
            ) = coefficient;
            directions(
                static_cast<Eigen::Index>(reaction_count + reaction),
                static_cast<Eigen::Index>(species_count + species)
            ) = coefficient;
        }
    }
    for (std::size_t species = 0; species < species_count; ++species) {
        directions(
            static_cast<Eigen::Index>(2 * reaction_count + species),
            static_cast<Eigen::Index>(species)
        ) = -1.0;
        directions(
            static_cast<Eigen::Index>(2 * reaction_count + species),
            static_cast<Eigen::Index>(species_count + species)
        ) = 1.0;
    }
    Eigen::FullPivLU<Eigen::MatrixXd> factor(directions);
    factor.setThreshold(
        256.0 * std::numeric_limits<double>::epsilon()
        * static_cast<double>(std::max(directions.rows(), directions.cols()))
    );
    return static_cast<std::size_t>(factor.rank());
}

struct PhaseBlock {
    std::vector<double> amounts;
    double log_packing = 0.0;
    double volume = 0.0;
    double packing = 0.0;
    double pressure = 0.0;
    double objective = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    std::vector<double> chemical_potentials;
    std::vector<double> physical_gradient;
    std::vector<double> physical_hessian;
    std::vector<double> state_jacobian;
    std::vector<double> objective_parameter_cross;
    std::vector<double> volume_parameter_direct;
    std::vector<double> pressure_parameter_direct;
    std::vector<double> chemical_potential_parameter_direct;
};

PhaseBlock evaluate_phase(
    const ProviderContext& provider,
    const ProviderActiveParameterSet& active_parameters,
    const std::vector<double>& g_ref,
    double temperature_k,
    double pressure_pa,
    const std::array<double, 2>& packing_interval,
    const std::vector<double>& variables,
    std::size_t offset
) {
    const std::size_t species_count = g_ref.size();
    const std::size_t phase_dimension = species_count + 1;
    const std::size_t parameter_count = active_parameters.identities.size();
    PhaseBlock result;
    result.amounts.resize(species_count);
    for (std::size_t species = 0; species < species_count; ++species) {
        result.amounts[species] = std::exp(variables[offset + species]);
        if (!std::isfinite(result.amounts[species]) || result.amounts[species] <= 0.0) {
            throw std::domain_error("reactive finite-phase amount chart is invalid");
        }
    }
    result.log_packing = variables[offset + species_count];
    result.packing = std::exp(result.log_packing);
    if (!std::isfinite(result.packing) || result.packing <= 0.0) {
        throw std::domain_error("reactive finite-phase packing chart is invalid");
    }

    const InversePackingGeometryEvaluation geometry =
        provider.evaluate_inverse_packing_geometry(
            temperature_k,
            result.amounts,
            result.log_packing,
            active_parameters
        );
    const std::size_t geometry_dimension = phase_dimension + parameter_count;
    if (geometry.gradient.size() != geometry_dimension
        || geometry.hessian.size() != geometry_dimension * geometry_dimension) {
        throw std::invalid_argument(
            "reactive finite-phase inverse-packing tensor dimensions changed"
        );
    }
    result.volume = geometry.volume_m3;

    MixturePhaseEvaluation mechanical;
    std::vector<double> state_parameter_cross;
    if (parameter_count == 0) {
        mechanical = provider.evaluate_mixture(
            temperature_k, result.amounts, result.volume
        );
    } else {
        const ParameterizedPhaseEvaluation parameterized =
            provider.evaluate_parameterized_phase(
                temperature_k,
                result.amounts,
                result.volume,
                packing_interval[0],
                packing_interval[1],
                active_parameters
            );
        mechanical = parameterized.phase;
        state_parameter_cross = parameterized.state_parameter_derivatives;
        result.pressure_parameter_direct =
            parameterized.pressure_parameter_derivatives_pa;
        result.chemical_potential_parameter_direct =
            parameterized.chemical_potential_parameter_derivatives_over_rt;
    }
    if (mechanical.gradient.size() != phase_dimension
        || mechanical.hessian.size() != phase_dimension * phase_dimension
        || (parameter_count != 0
            && state_parameter_cross.size() != phase_dimension * parameter_count)) {
        throw std::invalid_argument(
            "reactive finite-phase Provider tensor dimensions changed"
        );
    }
    result.pressure = mechanical.pressure_pa;
    result.physical_hessian = mechanical.hessian;
    result.physical_gradient = mechanical.gradient;
    const double pressure_over_rt = pressure_pa / (kGasConstantJPerMolK * temperature_k);
    for (std::size_t species = 0; species < species_count; ++species) {
        result.physical_gradient[species] += g_ref[species];
    }
    result.physical_gradient.back() += pressure_over_rt;
    result.objective = mechanical.value + pressure_over_rt * result.volume;
    for (std::size_t species = 0; species < species_count; ++species) {
        result.objective += g_ref[species] * result.amounts[species];
    }

    std::vector<double> physical_from_y(phase_dimension * phase_dimension, 0.0);
    for (std::size_t species = 0; species < species_count; ++species) {
        physical_from_y[species * phase_dimension + species] = 1.0;
    }
    for (std::size_t coordinate = 0; coordinate < phase_dimension; ++coordinate) {
        physical_from_y[species_count * phase_dimension + coordinate] =
            geometry.gradient[coordinate];
    }
    std::vector<double> y_gradient(phase_dimension, 0.0);
    std::vector<double> y_hessian(phase_dimension * phase_dimension, 0.0);
    for (std::size_t left = 0; left < phase_dimension; ++left) {
        for (std::size_t physical = 0; physical < phase_dimension; ++physical) {
            y_gradient[left] +=
                physical_from_y[physical * phase_dimension + left]
                * result.physical_gradient[physical];
        }
        for (std::size_t right = 0; right < phase_dimension; ++right) {
            for (std::size_t physical_left = 0;
                 physical_left < phase_dimension;
                 ++physical_left) {
                for (std::size_t physical_right = 0;
                     physical_right < phase_dimension;
                     ++physical_right) {
                    y_hessian[left * phase_dimension + right] +=
                        physical_from_y[physical_left * phase_dimension + left]
                        * mechanical.hessian[
                            physical_left * phase_dimension + physical_right
                        ]
                        * physical_from_y[
                            physical_right * phase_dimension + right
                        ];
                }
            }
            y_hessian[left * phase_dimension + right] +=
                result.physical_gradient.back()
                * geometry.hessian[left * geometry_dimension + right];
        }
    }

    std::vector<double> y_from_z(phase_dimension, 1.0);
    for (std::size_t species = 0; species < species_count; ++species) {
        y_from_z[species] = result.amounts[species];
    }
    result.gradient.resize(phase_dimension, 0.0);
    result.hessian.resize(phase_dimension * phase_dimension, 0.0);
    result.state_jacobian.assign(phase_dimension * phase_dimension, 0.0);
    for (std::size_t coordinate = 0; coordinate < phase_dimension; ++coordinate) {
        result.gradient[coordinate] = y_gradient[coordinate] * y_from_z[coordinate];
        for (std::size_t right = 0; right < phase_dimension; ++right) {
            result.hessian[coordinate * phase_dimension + right] =
                y_hessian[coordinate * phase_dimension + right]
                * y_from_z[coordinate] * y_from_z[right];
        }
        if (coordinate < species_count) {
            result.hessian[coordinate * phase_dimension + coordinate] +=
                y_gradient[coordinate] * result.amounts[coordinate];
        }
        for (std::size_t physical = 0; physical < phase_dimension; ++physical) {
            result.state_jacobian[physical * phase_dimension + coordinate] =
                physical_from_y[physical * phase_dimension + coordinate]
                * y_from_z[coordinate];
        }
    }

    result.chemical_potentials.resize(species_count, 0.0);
    for (std::size_t species = 0; species < species_count; ++species) {
        result.chemical_potentials[species] = mechanical.gradient[species] + g_ref[species];
    }

    result.objective_parameter_cross.assign(
        phase_dimension * parameter_count, 0.0
    );
    result.volume_parameter_direct.assign(parameter_count, 0.0);
    for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
        const std::size_t parameter_coordinate = phase_dimension + parameter;
        const double volume_direct = geometry.gradient[parameter_coordinate];
        result.volume_parameter_direct[parameter] = volume_direct;
        std::vector<double> physical_gradient_parameter(phase_dimension, 0.0);
        for (std::size_t physical = 0; physical < phase_dimension; ++physical) {
            physical_gradient_parameter[physical] =
                mechanical.hessian[physical * phase_dimension + species_count]
                    * volume_direct
                + state_parameter_cross[physical * parameter_count + parameter];
        }
        for (std::size_t y = 0; y < phase_dimension; ++y) {
            double derivative = geometry.hessian[
                y * geometry_dimension + parameter_coordinate
            ] * result.physical_gradient.back();
            for (std::size_t physical = 0; physical < phase_dimension; ++physical) {
                derivative +=
                    physical_from_y[physical * phase_dimension + y]
                    * physical_gradient_parameter[physical];
            }
            result.objective_parameter_cross[y * parameter_count + parameter] =
                y_from_z[y] * derivative;
        }
    }
    return result;
}

struct MultiphaseEvaluation {
    double objective = 0.0;
    std::vector<double> gradient;
    std::vector<double> constraints;
    std::vector<double> jacobian;
    std::vector<double> objective_hessian;
    std::vector<double> lagrangian_hessian;
    std::vector<PhaseBlock> phases;
};

MultiphaseEvaluation evaluate_multiphase(
    const ProviderContext& provider,
    const ProviderActiveParameterSet& active_parameters,
    const CompiledReactionSystem& system,
    const NormalizedBalances& balances,
    double temperature_k,
    double pressure_pa,
    const std::vector<std::array<double, 2>>& packing_intervals,
    const std::vector<double>& variables,
    const std::vector<double>& multipliers,
    std::size_t phase_count
) {
    const std::size_t species_count = system.species_count;
    const std::size_t phase_dimension = species_count + 1;
    const std::size_t variable_count = phase_count * phase_dimension;
    if (variables.size() != variable_count
        || multipliers.size() != balances.matrix.rows) {
        throw std::invalid_argument("reactive finite-phase NLP dimensions changed");
    }
    MultiphaseEvaluation result;
    result.gradient.assign(variable_count, 0.0);
    result.constraints.assign(balances.matrix.rows, 0.0);
    result.jacobian.assign(balances.matrix.rows * variable_count, 0.0);
    result.objective_hessian.assign(variable_count * variable_count, 0.0);
    result.phases.reserve(phase_count);
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        const std::size_t offset = phase * phase_dimension;
        result.phases.push_back(evaluate_phase(
            provider,
            active_parameters,
            system.g_ref,
            temperature_k,
            pressure_pa,
            packing_intervals[phase],
            variables,
            offset
        ));
        const PhaseBlock& block = result.phases.back();
        result.objective += block.objective;
        for (std::size_t coordinate = 0; coordinate < phase_dimension; ++coordinate) {
            result.gradient[offset + coordinate] = block.gradient[coordinate];
            for (std::size_t right = 0; right < phase_dimension; ++right) {
                result.objective_hessian[
                    (offset + coordinate) * variable_count + offset + right
                ] = block.hessian[coordinate * phase_dimension + right];
            }
        }
    }
    for (std::size_t row = 0; row < balances.matrix.rows; ++row) {
        result.constraints[row] = -balances.totals[row];
        for (std::size_t phase = 0; phase < phase_count; ++phase) {
            const std::size_t offset = phase * phase_dimension;
            const PhaseBlock& block = result.phases[phase];
            for (std::size_t species = 0; species < species_count; ++species) {
                result.constraints[row] +=
                    balances.matrix(row, species) * block.amounts[species];
                result.jacobian[row * variable_count + offset + species] =
                    balances.matrix(row, species) * block.amounts[species];
            }
        }
    }
    result.lagrangian_hessian = result.objective_hessian;
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        const std::size_t offset = phase * phase_dimension;
        for (std::size_t species = 0; species < species_count; ++species) {
            for (std::size_t row = 0; row < balances.matrix.rows; ++row) {
                result.lagrangian_hessian[
                    (offset + species) * variable_count + offset + species
                ] += multipliers[row]
                    * balances.matrix(row, species)
                    * result.phases[phase].amounts[species];
            }
        }
    }
    return result;
}

double relative_difference(double analytic, double numerical) {
    return std::abs(analytic - numerical)
        / std::max({1.0, std::abs(analytic), std::abs(numerical)});
}

std::vector<double> lagrangian_gradient(
    const MultiphaseEvaluation& evaluation,
    const std::vector<double>& multipliers
) {
    const std::size_t variable_count = evaluation.gradient.size();
    std::vector<double> result = evaluation.gradient;
    for (std::size_t row = 0; row < multipliers.size(); ++row) {
        for (std::size_t column = 0; column < variable_count; ++column) {
            result[column] += multipliers[row]
                * evaluation.jacobian[row * variable_count + column];
        }
    }
    return result;
}

bool audit_multiphase_derivatives(
    ReactiveFinitePhaseResult& result,
    const ProviderContext& provider,
    const ProviderActiveParameterSet& parameters,
    const CompiledReactionSystem& system,
    const NormalizedBalances& balances,
    double temperature_k,
    double pressure_pa,
    const std::vector<std::array<double, 2>>& packing_intervals,
    const std::vector<double>& solution,
    const std::vector<double>& multipliers,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    std::size_t phase_count
) {
    const std::size_t variable_count = solution.size();
    const MultiphaseEvaluation baseline = evaluate_multiphase(
        provider,
        parameters,
        system,
        balances,
        temperature_k,
        pressure_pa,
        packing_intervals,
        solution,
        multipliers,
        phase_count
    );
    double objective_error = 0.0;
    double constraint_error = 0.0;
    double hessian_error = 0.0;
    for (std::size_t column = 0; column < variable_count; ++column) {
        const double margin = std::min(
            solution[column] - lower[column], upper[column] - solution[column]
        );
        const double step = std::min(
            1.0e-5 * std::max(1.0, std::abs(solution[column])), 0.25 * margin
        );
        if (!std::isfinite(step) || step <= 1.0e-10) {
            result.derivative_audits_performed = true;
            result.derivative_audit_status = "failed";
            return false;
        }
        std::vector<double> plus_variables = solution;
        std::vector<double> minus_variables = solution;
        plus_variables[column] += step;
        minus_variables[column] -= step;
        const MultiphaseEvaluation plus = evaluate_multiphase(
            provider,
            parameters,
            system,
            balances,
            temperature_k,
            pressure_pa,
            packing_intervals,
            plus_variables,
            multipliers,
            phase_count
        );
        const MultiphaseEvaluation minus = evaluate_multiphase(
            provider,
            parameters,
            system,
            balances,
            temperature_k,
            pressure_pa,
            packing_intervals,
            minus_variables,
            multipliers,
            phase_count
        );
        const double inverse_span = 0.5 / step;
        objective_error = std::max(
            objective_error,
            relative_difference(
                baseline.gradient[column],
                (plus.objective - minus.objective) * inverse_span
            )
        );
        for (std::size_t row = 0; row < balances.matrix.rows; ++row) {
            constraint_error = std::max(
                constraint_error,
                relative_difference(
                    baseline.jacobian[row * variable_count + column],
                    (plus.constraints[row] - minus.constraints[row]) * inverse_span
                )
            );
        }
        const std::vector<double> plus_lagrangian_gradient =
            lagrangian_gradient(plus, multipliers);
        const std::vector<double> minus_lagrangian_gradient =
            lagrangian_gradient(minus, multipliers);
        for (std::size_t row = 0; row < variable_count; ++row) {
            hessian_error = std::max(
                hessian_error,
                relative_difference(
                    baseline.lagrangian_hessian[row * variable_count + column],
                    (plus_lagrangian_gradient[row] - minus_lagrangian_gradient[row])
                        * inverse_span
                )
            );
        }
    }
    result.derivative_audits_performed = true;
    result.objective_gradient_audit_max_relative_error = objective_error;
    result.constraint_jacobian_audit_max_relative_error = constraint_error;
    result.lagrangian_hessian_audit_max_relative_error = hessian_error;
    const bool passed = objective_error <= 2.0e-6
        && constraint_error <= 2.0e-8
        && hessian_error <= 2.0e-5;
    result.derivative_audit_status = passed ? "passed" : "failed";
    return passed;
}

class ReactiveFinitePhaseTnlp final : public Ipopt::TNLP {
public:
    ReactiveFinitePhaseTnlp(
        const ProviderContext& provider,
        ProviderActiveParameterSet active_parameters,
        const CompiledReactionSystem& system,
        NormalizedBalances balances,
        double temperature_k,
        double pressure_pa,
        std::vector<std::array<double, 2>> packing_intervals,
        std::vector<double> initial,
        std::vector<double> lower,
        std::vector<double> upper,
        std::size_t phase_count
    )
        : provider_(provider),
          active_parameters_(std::move(active_parameters)),
          system_(system),
          balances_(std::move(balances)),
          temperature_k_(temperature_k),
          pressure_pa_(pressure_pa),
          packing_intervals_(std::move(packing_intervals)),
          initial_(std::move(initial)),
          lower_(std::move(lower)),
          upper_(std::move(upper)),
          phase_count_(phase_count),
          solution_(initial_),
          multipliers_(balances_.matrix.rows, 0.0),
          lower_multipliers_(initial_.size(), 0.0),
          upper_multipliers_(initial_.size(), 0.0) {}

    bool get_nlp_info(
        Ipopt::Index& n,
        Ipopt::Index& m,
        Ipopt::Index& nnz_jac_g,
        Ipopt::Index& nnz_h_lag,
        IndexStyleEnum& index_style
    ) override {
        n = static_cast<Ipopt::Index>(initial_.size());
        m = static_cast<Ipopt::Index>(multipliers_.size());
        nnz_jac_g = n * m;
        nnz_h_lag = n * (n + 1) / 2;
        index_style = TNLP::C_STYLE;
        return true;
    }

    bool get_bounds_info(
        Ipopt::Index n,
        Ipopt::Number* x_l,
        Ipopt::Number* x_u,
        Ipopt::Index m,
        Ipopt::Number* g_l,
        Ipopt::Number* g_u
    ) override {
        if (!dimensions(n, m) || x_l == nullptr || x_u == nullptr
            || g_l == nullptr || g_u == nullptr) {
            return false;
        }
        std::copy(lower_.begin(), lower_.end(), x_l);
        std::copy(upper_.begin(), upper_.end(), x_u);
        std::fill(g_l, g_l + m, 0.0);
        std::fill(g_u, g_u + m, 0.0);
        return true;
    }

    bool get_starting_point(
        Ipopt::Index n,
        bool init_x,
        Ipopt::Number* x,
        bool,
        Ipopt::Number*,
        Ipopt::Number*,
        Ipopt::Index m,
        bool,
        Ipopt::Number*
    ) override {
        if (!dimensions(n, m) || !init_x || x == nullptr) {
            return false;
        }
        std::copy(initial_.begin(), initial_.end(), x);
        return true;
    }

    bool eval_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number& value
    ) override {
        if (n != static_cast<Ipopt::Index>(initial_.size()) || x == nullptr) {
            return false;
        }
        return guarded([&] {
            value = evaluate(x, zero_multipliers()).objective;
        });
    }

    bool eval_grad_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number* gradient
    ) override {
        if (n != static_cast<Ipopt::Index>(initial_.size())
            || x == nullptr || gradient == nullptr) {
            return false;
        }
        return guarded([&] {
            const MultiphaseEvaluation evaluated = evaluate(x, zero_multipliers());
            std::copy(evaluated.gradient.begin(), evaluated.gradient.end(), gradient);
        });
    }

    bool eval_g(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Index m,
        Ipopt::Number* constraints
    ) override {
        if (!dimensions(n, m) || x == nullptr || constraints == nullptr) {
            return false;
        }
        return guarded([&] {
            const MultiphaseEvaluation evaluated = evaluate(x, zero_multipliers());
            std::copy(
                evaluated.constraints.begin(), evaluated.constraints.end(), constraints
            );
        });
    }

    bool eval_jac_g(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Index m,
        Ipopt::Index nonzero_count,
        Ipopt::Index* rows,
        Ipopt::Index* columns,
        Ipopt::Number* values
    ) override {
        if (!dimensions(n, m) || nonzero_count != n * m) {
            return false;
        }
        if (values == nullptr) {
            if (rows == nullptr || columns == nullptr) {
                return false;
            }
            for (Ipopt::Index row = 0; row < m; ++row) {
                for (Ipopt::Index column = 0; column < n; ++column) {
                    rows[row * n + column] = row;
                    columns[row * n + column] = column;
                }
            }
            return true;
        }
        if (x == nullptr) {
            return false;
        }
        return guarded([&] {
            const MultiphaseEvaluation evaluated = evaluate(x, zero_multipliers());
            std::copy(evaluated.jacobian.begin(), evaluated.jacobian.end(), values);
        });
    }

    bool eval_h(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number objective_factor,
        Ipopt::Index m,
        const Ipopt::Number* multipliers,
        bool,
        Ipopt::Index nonzero_count,
        Ipopt::Index* rows,
        Ipopt::Index* columns,
        Ipopt::Number* values
    ) override {
        if (!dimensions(n, m) || nonzero_count != n * (n + 1) / 2) {
            return false;
        }
        if (values == nullptr) {
            if (rows == nullptr || columns == nullptr) {
                return false;
            }
            Ipopt::Index entry = 0;
            for (Ipopt::Index row = 0; row < n; ++row) {
                for (Ipopt::Index column = 0; column <= row; ++column) {
                    rows[entry] = row;
                    columns[entry] = column;
                    ++entry;
                }
            }
            return true;
        }
        if (x == nullptr || multipliers == nullptr) {
            return false;
        }
        return guarded([&] {
            const std::vector<double> lambda(multipliers, multipliers + m);
            const MultiphaseEvaluation evaluated = evaluate(x, lambda);
            Ipopt::Index entry = 0;
            for (Ipopt::Index row = 0; row < n; ++row) {
                for (Ipopt::Index column = 0; column <= row; ++column) {
                    const std::size_t index = static_cast<std::size_t>(row * n + column);
                    values[entry++] = objective_factor
                            * evaluated.objective_hessian[index]
                        + evaluated.lagrangian_hessian[index]
                        - evaluated.objective_hessian[index];
                }
            }
        });
    }

    bool intermediate_callback(
        Ipopt::AlgorithmMode,
        Ipopt::Index iteration,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Index,
        const Ipopt::IpoptData*,
        Ipopt::IpoptCalculatedQuantities*
    ) override {
        iterations_ = std::max(
            iterations_, static_cast<std::size_t>(std::max<Ipopt::Index>(0, iteration))
        );
        return callback_error_.empty();
    }

    void finalize_solution(
        Ipopt::SolverReturn,
        Ipopt::Index n,
        const Ipopt::Number* x,
        const Ipopt::Number* z_l,
        const Ipopt::Number* z_u,
        Ipopt::Index m,
        const Ipopt::Number*,
        const Ipopt::Number* lambda,
        Ipopt::Number,
        const Ipopt::IpoptData*,
        Ipopt::IpoptCalculatedQuantities*
    ) override {
        if (n == static_cast<Ipopt::Index>(solution_.size())
            && x != nullptr && z_l != nullptr && z_u != nullptr) {
            std::copy(x, x + n, solution_.begin());
            std::copy(z_l, z_l + n, lower_multipliers_.begin());
            std::copy(z_u, z_u + n, upper_multipliers_.begin());
        }
        if (m == static_cast<Ipopt::Index>(multipliers_.size()) && lambda != nullptr) {
            std::copy(lambda, lambda + m, multipliers_.begin());
        }
    }

    [[nodiscard]] const std::vector<double>& solution() const { return solution_; }
    [[nodiscard]] const std::vector<double>& multipliers() const { return multipliers_; }
    [[nodiscard]] const std::vector<double>& lower_multipliers() const {
        return lower_multipliers_;
    }
    [[nodiscard]] const std::vector<double>& upper_multipliers() const {
        return upper_multipliers_;
    }
    [[nodiscard]] const std::string& callback_error() const { return callback_error_; }
    [[nodiscard]] const std::string& callback_failure_kind() const {
        return callback_failure_kind_;
    }
    [[nodiscard]] std::size_t iterations() const { return iterations_; }

private:
    [[nodiscard]] bool dimensions(Ipopt::Index n, Ipopt::Index m) const {
        return n == static_cast<Ipopt::Index>(initial_.size())
            && m == static_cast<Ipopt::Index>(multipliers_.size());
    }

    [[nodiscard]] std::vector<double> zero_multipliers() const {
        return std::vector<double>(multipliers_.size(), 0.0);
    }

    MultiphaseEvaluation evaluate(
        const Ipopt::Number* x,
        const std::vector<double>& multipliers
    ) const {
        return evaluate_multiphase(
            provider_,
            active_parameters_,
            system_,
            balances_,
            temperature_k_,
            pressure_pa_,
            packing_intervals_,
            std::vector<double>(x, x + initial_.size()),
            multipliers,
            phase_count_
        );
    }

    template <class Callback>
    bool guarded(Callback&& callback) {
        try {
            callback();
            return true;
        } catch (const ProviderDomainError& error) {
            callback_error_ = error.what();
            callback_failure_kind_ = "provider_domain";
            return false;
        } catch (const std::domain_error& error) {
            callback_error_ = error.what();
            callback_failure_kind_ = "physical_domain";
            return false;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            callback_failure_kind_ = "internal_contract";
            return false;
        }
    }

    const ProviderContext& provider_;
    ProviderActiveParameterSet active_parameters_;
    const CompiledReactionSystem& system_;
    NormalizedBalances balances_;
    double temperature_k_;
    double pressure_pa_;
    std::vector<std::array<double, 2>> packing_intervals_;
    std::vector<double> initial_;
    std::vector<double> lower_;
    std::vector<double> upper_;
    std::size_t phase_count_;
    std::vector<double> solution_;
    std::vector<double> multipliers_;
    std::vector<double> lower_multipliers_;
    std::vector<double> upper_multipliers_;
    std::string callback_error_;
    std::string callback_failure_kind_ = "none";
    std::size_t iterations_ = 0;
};

double phase_distinctness(const std::vector<PhaseBlock>& phases) {
    double result = 0.0;
    for (std::size_t left = 0; left < phases.size(); ++left) {
        const double left_total = std::accumulate(
            phases[left].amounts.begin(), phases[left].amounts.end(), 0.0
        );
        for (std::size_t right = left + 1; right < phases.size(); ++right) {
            const double right_total = std::accumulate(
                phases[right].amounts.begin(), phases[right].amounts.end(), 0.0
            );
            double distance = 0.0;
            for (std::size_t species = 0; species < phases[left].amounts.size(); ++species) {
                distance = std::max(
                    distance,
                    std::abs(std::log(
                        (phases[left].amounts[species] / left_total)
                        / (phases[right].amounts[species] / right_total)
                    ))
                );
            }
            result = std::max(result, distance);
        }
    }
    return result;
}

double reduced_hessian_minimum(
    const MultiphaseEvaluation& evaluation,
    std::size_t variable_count,
    std::size_t constraint_count
) {
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
        jacobian(
            evaluation.jacobian.data(),
            static_cast<Eigen::Index>(constraint_count),
            static_cast<Eigen::Index>(variable_count)
        );
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
        hessian(
            evaluation.lagrangian_hessian.data(),
            static_cast<Eigen::Index>(variable_count),
            static_cast<Eigen::Index>(variable_count)
        );
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(jacobian, Eigen::ComputeFullV);
    const double threshold = 256.0 * std::numeric_limits<double>::epsilon()
        * static_cast<double>(std::max(variable_count, constraint_count))
        * (svd.singularValues().size() == 0 ? 1.0 : svd.singularValues()(0));
    std::size_t rank = 0;
    for (Eigen::Index index = 0; index < svd.singularValues().size(); ++index) {
        rank += static_cast<std::size_t>(svd.singularValues()(index) > threshold);
    }
    if (rank >= variable_count) {
        return std::numeric_limits<double>::infinity();
    }
    const Eigen::MatrixXd nullspace = svd.matrixV().rightCols(
        static_cast<Eigen::Index>(variable_count - rank)
    );
    const Eigen::MatrixXd reduced = nullspace.transpose()
        * (0.5 * (hessian + hessian.transpose())) * nullspace;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(reduced);
    if (eigen.info() != Eigen::Success) {
        throw std::runtime_error("reactive finite-phase reduced Hessian eigensolve failed");
    }
    return eigen.eigenvalues().minCoeff();
}

void populate_sensitivities(
    ReactiveFinitePhaseResult& result,
    const MultiphaseEvaluation& evaluation,
    const ProviderActiveParameterSet& active_parameters,
    const std::vector<double>& multipliers,
    double temperature_k,
    std::size_t species_count
) {
    const std::size_t parameter_count = active_parameters.identities.size();
    if (parameter_count == 0) {
        result.sensitivity_status = "not_requested";
        return;
    }
    const std::size_t phase_count = evaluation.phases.size();
    const std::size_t phase_dimension = species_count + 1;
    const std::size_t variable_count = phase_count * phase_dimension;
    const std::size_t constraint_count = evaluation.constraints.size();
    const std::size_t kkt_dimension = variable_count + constraint_count;
    Eigen::MatrixXd kkt = Eigen::MatrixXd::Zero(
        static_cast<Eigen::Index>(kkt_dimension),
        static_cast<Eigen::Index>(kkt_dimension)
    );
    for (std::size_t row = 0; row < variable_count; ++row) {
        for (std::size_t column = 0; column < variable_count; ++column) {
            kkt(
                static_cast<Eigen::Index>(row),
                static_cast<Eigen::Index>(column)
            ) = evaluation.lagrangian_hessian[row * variable_count + column];
        }
    }
    for (std::size_t row = 0; row < constraint_count; ++row) {
        for (std::size_t column = 0; column < variable_count; ++column) {
            const double value = evaluation.jacobian[row * variable_count + column];
            kkt(
                static_cast<Eigen::Index>(column),
                static_cast<Eigen::Index>(variable_count + row)
            ) = value;
            kkt(
                static_cast<Eigen::Index>(variable_count + row),
                static_cast<Eigen::Index>(column)
            ) = value;
        }
    }
    Eigen::FullPivLU<Eigen::MatrixXd> factor(kkt);
    if (static_cast<std::size_t>(factor.rank()) != kkt_dimension) {
        result.sensitivity_status = "unavailable";
        result.sensitivity_failure_reason = "singular_exact_kkt_system";
        return;
    }
    result.parameter_order = active_parameters.identities;
    result.amount_derivatives.assign(
        phase_count * species_count * parameter_count, 0.0
    );
    result.volume_derivatives.assign(phase_count * parameter_count, 0.0);
    result.pressure_derivatives.assign(phase_count * parameter_count, 0.0);
    result.chemical_potential_derivatives.assign(
        phase_count * species_count * parameter_count, 0.0
    );
    for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
        Eigen::VectorXd right = Eigen::VectorXd::Zero(
            static_cast<Eigen::Index>(kkt_dimension)
        );
        for (std::size_t phase = 0; phase < phase_count; ++phase) {
            const std::size_t offset = phase * phase_dimension;
            for (std::size_t coordinate = 0; coordinate < phase_dimension; ++coordinate) {
                right(static_cast<Eigen::Index>(offset + coordinate)) =
                    -evaluation.phases[phase].objective_parameter_cross[
                        coordinate * parameter_count + parameter
                    ];
            }
        }
        const Eigen::VectorXd response = factor.solve(right);
        if (!response.allFinite()) {
            result.sensitivity_status = "unavailable";
            result.sensitivity_failure_reason = "nonfinite_exact_kkt_response";
            return;
        }
        for (std::size_t phase = 0; phase < phase_count; ++phase) {
            const PhaseBlock& block = evaluation.phases[phase];
            const std::size_t offset = phase * phase_dimension;
            std::vector<double> physical_response(phase_dimension, 0.0);
            for (std::size_t physical = 0; physical < phase_dimension; ++physical) {
                for (std::size_t coordinate = 0; coordinate < phase_dimension; ++coordinate) {
                    physical_response[physical] +=
                        block.state_jacobian[physical * phase_dimension + coordinate]
                        * response(static_cast<Eigen::Index>(offset + coordinate));
                }
            }
            physical_response.back() += block.volume_parameter_direct[parameter];
            for (std::size_t species = 0; species < species_count; ++species) {
                result.amount_derivatives[
                    (phase * species_count + species) * parameter_count + parameter
                ] = physical_response[species];
                double chemical_response =
                    block.chemical_potential_parameter_direct[
                        species * parameter_count + parameter
                    ];
                for (std::size_t physical = 0; physical < phase_dimension; ++physical) {
                    chemical_response +=
                        block.physical_hessian[species * phase_dimension + physical]
                        * physical_response[physical];
                }
                result.chemical_potential_derivatives[
                    (phase * species_count + species) * parameter_count + parameter
                ] = chemical_response;
            }
            result.volume_derivatives[phase * parameter_count + parameter] =
                physical_response.back();
            double pressure_response = block.pressure_parameter_direct[parameter];
            for (std::size_t physical = 0; physical < phase_dimension; ++physical) {
                pressure_response -= kGasConstantJPerMolK * temperature_k
                    * block.physical_hessian[
                        species_count * phase_dimension + physical
                    ] * physical_response[physical];
            }
            result.pressure_derivatives[phase * parameter_count + parameter] =
                pressure_response;
        }
    }
    result.sensitivity_status = "available";
    result.sensitivity_failure_reason.clear();
    (void)multipliers;
}

}  // namespace

ReactiveFinitePhaseResult solve_provider_reactive_finite_phases(
    const CompiledReactionSystem& system,
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    const std::vector<std::array<double, 2>>& packing_fraction_intervals,
    const std::vector<ReactiveFinitePhaseStart>& starts,
    double trace_floor,
    const ProviderActiveParameterSet* active_parameters,
    int max_iterations,
    bool perform_derivative_audits
) {
    ReactiveFinitePhaseResult result;
    const std::size_t phase_count = starts.size();
    const std::size_t species_count = system.species_count;
    const std::size_t phase_dimension = species_count + 1;
    const std::size_t variable_count = phase_count * phase_dimension;
    if (phase_count != 2 || packing_fraction_intervals.size() != phase_count) {
        throw std::invalid_argument(
            "reactive finite-phase admission requires exactly two declared phases"
        );
    }
    if (!std::all_of(system.charges.begin(), system.charges.end(), [](int charge) {
            return charge == 0;
        })) {
        throw std::invalid_argument(
            "reactive finite-phase admission currently requires neutral species"
        );
    }
    if (system.species_count != system.original_species_count
        || !system.removed_species_indices.empty()) {
        throw std::invalid_argument(
            "reactive finite-phase admission requires a complete strict-interior species face"
        );
    }
    if (!std::isfinite(temperature_k) || temperature_k <= 0.0
        || !std::isfinite(pressure_pa) || pressure_pa <= 0.0
        || !std::isfinite(trace_floor) || trace_floor <= 0.0
        || max_iterations <= 0) {
        throw std::invalid_argument("reactive finite-phase state controls are invalid");
    }
    result.admissible_nullity = 2 * species_count - system.balance_matrix.rows;
    result.declared_direction_rank = declared_direction_rank(system);
    if (result.declared_direction_rank != result.admissible_nullity) {
        throw std::invalid_argument(
            "declared reaction and transfer directions do not span the conservation nullspace"
        );
    }

    const ProviderActiveParameterSet parameters = active_parameters == nullptr
        ? ProviderActiveParameterSet{}
        : *active_parameters;
    std::vector<double> initial(variable_count, 0.0);
    std::vector<double> lower(variable_count, 0.0);
    std::vector<double> upper(variable_count, 0.0);
    // Positive molar masses make this a rigorous per-species amount ceiling.
    const double total_feed_mass_kg = std::inner_product(
        system.feed_amounts.begin(),
        system.feed_amounts.end(),
        system.molar_masses_kg_per_mol.begin(),
        0.0
    );
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        const ReactiveFinitePhaseStart& start = starts[phase];
        const std::array<double, 2>& packing = packing_fraction_intervals[phase];
        if (start.amounts_mol.size() != species_count
            || !std::all_of(start.amounts_mol.begin(), start.amounts_mol.end(),
                [trace_floor](double amount) {
                    return std::isfinite(amount) && amount > trace_floor;
                })
            || !std::isfinite(start.molar_volume_m3_per_mol)
            || start.molar_volume_m3_per_mol <= 0.0
            || !std::isfinite(packing[0]) || !std::isfinite(packing[1])
            || packing[0] <= 0.0 || packing[1] <= packing[0] || packing[1] >= 1.0) {
            throw std::invalid_argument("reactive finite-phase start or domain is invalid");
        }
        const std::size_t offset = phase * phase_dimension;
        const double total = std::accumulate(
            start.amounts_mol.begin(), start.amounts_mol.end(), 0.0
        );
        const double volume = total * start.molar_volume_m3_per_mol;
        const double start_packing = provider.evaluate_packing_fraction(
            temperature_k, start.amounts_mol, volume
        ).value;
        if (start_packing <= packing[0] || start_packing >= packing[1]) {
            throw std::invalid_argument(
                "reactive finite-phase start lies outside the declared packing interval"
            );
        }
        for (std::size_t species = 0; species < species_count; ++species) {
            initial[offset + species] = std::log(start.amounts_mol[species]);
            lower[offset + species] = std::log(0.1 * trace_floor);
            upper[offset + species] = std::log(std::nextafter(
                total_feed_mass_kg / system.molar_masses_kg_per_mol[species],
                std::numeric_limits<double>::infinity()
            ));
        }
        initial[offset + species_count] = std::log(start_packing);
        lower[offset + species_count] = std::log(packing[0]);
        upper[offset + species_count] = std::nextafter(
            std::log(packing[1]), -std::numeric_limits<double>::infinity()
        );
    }

    const NormalizedBalances balances = normalized_balances(system);
    const MultiphaseEvaluation initial_evaluation = evaluate_multiphase(
        provider,
        parameters,
        system,
        balances,
        temperature_k,
        pressure_pa,
        packing_fraction_intervals,
        initial,
        std::vector<double>(balances.matrix.rows, 0.0),
        phase_count
    );
    if (inf_norm(initial_evaluation.constraints) > 1.0e-8) {
        throw std::invalid_argument(
            "reactive finite-phase starts do not preserve the global conserved totals"
        );
    }

    ReactiveFinitePhaseTnlp* raw_problem = new ReactiveFinitePhaseTnlp(
        provider,
        parameters,
        system,
        balances,
        temperature_k,
        pressure_pa,
        packing_fraction_intervals,
        initial,
        lower,
        upper,
        phase_count
    );
    Ipopt::SmartPtr<Ipopt::TNLP> problem = raw_problem;
    Ipopt::SmartPtr<Ipopt::IpoptApplication> application = IpoptApplicationFactory();
    application->Options()->SetIntegerValue("print_level", 0);
    application->Options()->SetStringValue("sb", "yes");
    application->Options()->SetStringValue("hessian_approximation", "exact");
    application->Options()->SetStringValue("mu_strategy", "adaptive");
    application->Options()->SetNumericValue("tol", 1.0e-10);
    application->Options()->SetNumericValue("constr_viol_tol", 1.0e-10);
    application->Options()->SetNumericValue("acceptable_tol", 1.0e-8);
    application->Options()->SetIntegerValue("max_iter", max_iterations);
    const Ipopt::ApplicationReturnStatus initialized = application->Initialize();
    if (initialized != Ipopt::Solve_Succeeded) {
        result.solver_status = "initialization_" + ipopt_status_name(initialized);
        result.failure_reason = "ipopt_initialization_failed";
        result.numerical_status = "failed";
        return result;
    }
    const Ipopt::ApplicationReturnStatus status = application->OptimizeTNLP(problem);
    result.solver_status = ipopt_status_name(status);
    result.callback_error = raw_problem->callback_error();
    result.callback_failure_kind = raw_problem->callback_failure_kind();
    result.optimizer_iterations = raw_problem->iterations();
    result.provider_evaluation_count = static_cast<std::size_t>(provider.evaluation_count());
    if (status != Ipopt::Solve_Succeeded || !result.callback_error.empty()) {
        if (!result.callback_error.empty()) {
            if (result.callback_failure_kind == "provider_domain") {
                result.failure_reason = "provider_callback_domain_failed";
                result.provider_domain_status = "failed";
            } else {
                result.failure_reason = result.callback_failure_kind == "physical_domain"
                    ? "physical_domain_callback_failed"
                    : "internal_callback_failed";
                result.provider_domain_status = result.callback_failure_kind == "physical_domain"
                    ? "passed"
                    : "not_adjudicated";
                if (result.callback_failure_kind == "physical_domain") {
                    result.physical_status = "failed";
                } else {
                    result.numerical_status = "failed";
                }
            }
        } else {
            result.failure_reason = "ipopt_solve_failed";
            result.provider_domain_status = "passed";
            result.numerical_status = "failed";
        }
        return result;
    }
    result.provider_domain_status = "passed";
    result.numerical_status = "passed";

    const std::vector<double>& solution = raw_problem->solution();
    const std::vector<double>& multipliers = raw_problem->multipliers();
    const MultiphaseEvaluation evaluated = evaluate_multiphase(
        provider,
        parameters,
        system,
        balances,
        temperature_k,
        pressure_pa,
        packing_fraction_intervals,
        solution,
        multipliers,
        phase_count
    );
    result.objective_over_rt = evaluated.objective;
    result.balance_inf_norm = inf_norm(evaluated.constraints);
    result.phase_distinctness = phase_distinctness(evaluated.phases);
    result.phase_amounts_mol.reserve(phase_count);
    result.phase_volumes_m3.reserve(phase_count);
    result.phase_pressures_pa.reserve(phase_count);
    result.phase_packing_fractions.reserve(phase_count);
    result.phase_chemical_potentials_over_rt.reserve(phase_count);
    std::vector<double> stationarity = evaluated.gradient;
    for (std::size_t row = 0; row < balances.matrix.rows; ++row) {
        for (std::size_t column = 0; column < variable_count; ++column) {
            stationarity[column] += multipliers[row]
                * evaluated.jacobian[row * variable_count + column];
        }
    }
    result.kkt_stationarity_inf_norm = inf_norm(stationarity);
    result.pressure_relative_inf_norm = 0.0;
    result.reaction_affinity_inf_norm = 0.0;
    result.transfer_chemical_potential_inf_norm = 0.0;
    double minimum_amount = std::numeric_limits<double>::infinity();
    double minimum_bound_margin = std::numeric_limits<double>::infinity();
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        const PhaseBlock& block = evaluated.phases[phase];
        result.phase_amounts_mol.push_back(block.amounts);
        result.phase_volumes_m3.push_back(block.volume);
        result.phase_pressures_pa.push_back(block.pressure);
        result.phase_packing_fractions.push_back(block.packing);
        result.phase_chemical_potentials_over_rt.push_back(block.chemical_potentials);
        result.pressure_relative_inf_norm = std::max(
            result.pressure_relative_inf_norm,
            std::abs(block.pressure - pressure_pa) / pressure_pa
        );
        for (double amount : block.amounts) {
            minimum_amount = std::min(minimum_amount, amount);
        }
        const std::size_t offset = phase * phase_dimension;
        for (std::size_t coordinate = 0; coordinate < phase_dimension; ++coordinate) {
            minimum_bound_margin = std::min({
                minimum_bound_margin,
                solution[offset + coordinate] - lower[offset + coordinate],
                upper[offset + coordinate] - solution[offset + coordinate],
            });
        }
        for (std::size_t reaction = 0;
             reaction < system.reaction_matrix.rows;
             ++reaction) {
            double affinity = 0.0;
            for (std::size_t species = 0; species < species_count; ++species) {
                affinity += system.reaction_matrix(reaction, species)
                    * block.chemical_potentials[species];
            }
            result.reaction_affinity_inf_norm = std::max(
                result.reaction_affinity_inf_norm, std::abs(affinity)
            );
        }
    }
    for (std::size_t species = 0; species < species_count; ++species) {
        result.transfer_chemical_potential_inf_norm = std::max(
            result.transfer_chemical_potential_inf_norm,
            std::abs(
                evaluated.phases[0].chemical_potentials[species]
                - evaluated.phases[1].chemical_potentials[species]
            )
        );
    }
    result.reduced_hessian_min_eigenvalue = reduced_hessian_minimum(
        evaluated, variable_count, balances.matrix.rows
    );
    const double hessian_scale = std::max(
        1.0, inf_norm(evaluated.lagrangian_hessian)
    );
    const double strict_minimum_threshold =
        1.0e-10 * hessian_scale;
    const bool physical_certificates_passed =
        result.balance_inf_norm <= kChemicalBalanceTolerance
        && result.pressure_relative_inf_norm <= kChemicalPressureTolerance
        && result.reaction_affinity_inf_norm <= kChemicalAffinityTolerance
        && result.transfer_chemical_potential_inf_norm <= kChemicalAffinityTolerance
        && result.kkt_stationarity_inf_norm <= kChemicalKktTolerance
        && result.phase_distinctness > 1.0e-6
        && minimum_amount > trace_floor
        && minimum_bound_margin > 1.0e-7
        && result.reduced_hessian_min_eigenvalue > strict_minimum_threshold;
    result.physical_status = physical_certificates_passed ? "passed" : "failed";
    if (!physical_certificates_passed) {
        result.failure_reason = "reactive_finite_phase_certificate_failed";
        return result;
    }
    if (perform_derivative_audits) {
        bool derivative_audits_passed = false;
        try {
            derivative_audits_passed = audit_multiphase_derivatives(
                result,
                provider,
                parameters,
                system,
                balances,
                temperature_k,
                pressure_pa,
                packing_fraction_intervals,
                solution,
                multipliers,
                lower,
                upper,
                phase_count
            );
        } catch (const std::exception& error) {
            result.derivative_audits_performed = true;
            result.derivative_audit_status = "failed";
            result.callback_error = error.what();
        }
        if (!derivative_audits_passed) {
            result.numerical_status = "failed";
            result.failure_reason = "reactive_finite_phase_derivative_audit_failed";
            result.provider_evaluation_count =
                static_cast<std::size_t>(provider.evaluation_count());
            return result;
        }
    }
    result.accepted = true;
    populate_sensitivities(
        result,
        evaluated,
        parameters,
        multipliers,
        temperature_k,
        species_count
    );
    if (!parameters.identities.empty() && result.sensitivity_status != "available") {
        result.accepted = false;
        result.failure_reason = "reactive_finite_phase_sensitivity_failed";
    }
    result.provider_evaluation_count = static_cast<std::size_t>(provider.evaluation_count());
    return result;
}

}  // namespace epcsaft_engine::equilibrium
