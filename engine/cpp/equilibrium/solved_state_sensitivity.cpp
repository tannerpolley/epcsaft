#include "solved_state_sensitivity.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

namespace epcsaft_engine::equilibrium {
namespace {

constexpr double kGasConstantJPerMolK = 8.31446261815324;
constexpr double kPackingFractionMinimum = 1.0e-6;
constexpr double kPackingFractionMaximum = 0.74;
constexpr double kChemicalPotentialResidualMax = 1.0e-6;
constexpr double kPressureResidualMax = 1.0e-8;

bool finite(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
    });
}

std::size_t amount_state_index(
    std::size_t phase,
    std::size_t component,
    std::size_t component_count
) {
    return phase * component_count + component;
}

std::size_t volume_state_index(
    std::size_t phase,
    std::size_t phase_count,
    std::size_t component_count
) {
    return (phase_count - 1) * component_count + phase;
}

std::size_t output_amount_index(
    std::size_t phase,
    std::size_t component,
    std::size_t component_count
) {
    return phase * (component_count + 1) + component;
}

std::size_t output_volume_index(
    std::size_t phase,
    std::size_t component_count
) {
    return phase * (component_count + 1) + component_count;
}

void append_identities(
    SolvedStateSensitivity& result,
    const std::vector<std::string>& component_ids,
    std::size_t phase_count
) {
    const std::size_t reference = phase_count - 1;
    for (std::size_t phase = 0; phase < reference; ++phase) {
        for (const auto& component : component_ids) {
            result.state_identities.push_back(
                "phase[" + std::to_string(phase) + "].amount[" + component + "]"
            );
            result.residual_identities.push_back(
                "phase[" + std::to_string(phase) + "].mu[" + component
                + "]-phase[" + std::to_string(reference) + "].mu["
                + component + "]"
            );
        }
    }
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        result.state_identities.push_back(
            "phase[" + std::to_string(phase) + "].volume_m3"
        );
        result.residual_identities.push_back(
            "phase[" + std::to_string(phase) + "].pressure_over_target_minus_one"
        );
        for (const auto& component : component_ids) {
            result.output_identities.push_back(
                "phase[" + std::to_string(phase) + "].amount[" + component + "]"
            );
        }
        result.output_identities.push_back(
            "phase[" + std::to_string(phase) + "].volume_m3"
        );
    }
}

bool validate_phases(
    const std::vector<double>& feed,
    const std::vector<SolvedPhaseState>& phases,
    SolvedStateSensitivity& result
) {
    const std::size_t component_count = feed.size();
    if (component_count == 0 || phases.empty() || phases.size() > 3
        || !finite(feed)) {
        result.failure_reason = "unsupported_solved_topology";
        return false;
    }
    std::vector<double> reconstructed(component_count, 0.0);
    for (const auto& phase : phases) {
        if (phase.amounts_mol.size() != component_count
            || !finite(phase.amounts_mol)
            || !std::isfinite(phase.volume_m3) || phase.volume_m3 <= 0.0) {
            result.failure_reason = "invalid_solved_phase";
            return false;
        }
        const double total = std::accumulate(
            phase.amounts_mol.begin(), phase.amounts_mol.end(), 0.0
        );
        const double boundary = 4096.0 * std::numeric_limits<double>::epsilon()
            * std::max(1.0, std::abs(total));
        if (!std::isfinite(total) || total <= boundary
            || std::any_of(
                phase.amounts_mol.begin(),
                phase.amounts_mol.end(),
                [boundary](double amount) { return amount <= boundary; }
            )) {
            result.failure_reason = "phase_boundary_not_differentiable";
            return false;
        }
        for (std::size_t component = 0; component < component_count; ++component) {
            reconstructed[component] += phase.amounts_mol[component];
        }
    }
    for (std::size_t component = 0; component < component_count; ++component) {
        const double scale = std::max(1.0, std::abs(feed[component]));
        if (feed[component] <= 0.0
            || std::abs(reconstructed[component] - feed[component])
                > 2.0e-10 * scale) {
            result.failure_reason = "material_balance_mismatch";
            return false;
        }
    }
    for (std::size_t left = 0; left < phases.size(); ++left) {
        const double left_total = std::accumulate(
            phases[left].amounts_mol.begin(), phases[left].amounts_mol.end(), 0.0
        );
        for (std::size_t right = left + 1; right < phases.size(); ++right) {
            const double right_total = std::accumulate(
                phases[right].amounts_mol.begin(), phases[right].amounts_mol.end(), 0.0
            );
            double composition_distance = 0.0;
            for (std::size_t component = 0; component < component_count; ++component) {
                composition_distance = std::max(
                    composition_distance,
                    std::abs(
                        phases[left].amounts_mol[component] / left_total
                        - phases[right].amounts_mol[component] / right_total
                    )
                );
            }
            const double density_left = left_total / phases[left].volume_m3;
            const double density_right = right_total / phases[right].volume_m3;
            const double density_distance = std::abs(density_left - density_right)
                / std::max({1.0, std::abs(density_left), std::abs(density_right)});
            if (composition_distance <= 1.0e-9 && density_distance <= 1.0e-9) {
                result.failure_reason = "nonunique_coalesced_phases";
                return false;
            }
        }
    }
    return true;
}

bool validate_defining_equations(
    const std::vector<ParameterizedPhaseEvaluation>& evaluations,
    double pressure_pa,
    std::size_t component_count,
    SolvedStateSensitivity& result
) {
    const std::size_t reference = evaluations.size() - 1;
    for (std::size_t phase = 0; phase < reference; ++phase) {
        for (std::size_t component = 0; component < component_count; ++component) {
            const double residual = evaluations[phase].phase.gradient[component]
                - evaluations[reference].phase.gradient[component];
            if (!std::isfinite(residual)
                || std::abs(residual) > kChemicalPotentialResidualMax) {
                result.failure_reason = "defining_equations_not_satisfied";
                return false;
            }
        }
    }
    for (const auto& evaluation : evaluations) {
        const double residual = evaluation.phase.pressure_pa / pressure_pa - 1.0;
        if (!std::isfinite(residual)
            || std::abs(residual) > kPressureResidualMax) {
            result.failure_reason = "defining_equations_not_satisfied";
            return false;
        }
    }
    return true;
}

struct RankRevealingSolve {
    std::size_t rank = 0;
    double condition_number = std::numeric_limits<double>::infinity();
    double backward_error = std::numeric_limits<double>::infinity();
    std::vector<double> solution;
};

RankRevealingSolve solve_total_derivatives(
    const std::vector<double>& f_y,
    const std::vector<double>& f_p,
    std::size_t dimension,
    std::size_t parameter_count
) {
    RankRevealingSolve result;
    if (dimension == 0 || parameter_count == 0
        || f_y.size() != dimension * dimension
        || f_p.size() != dimension * parameter_count
        || !finite(f_y) || !finite(f_p)) {
        return result;
    }
    using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    Eigen::Map<const Matrix> matrix(f_y.data(), dimension, dimension);
    Eigen::Map<const Matrix> partials(f_p.data(), dimension, parameter_count);
    Eigen::VectorXd row_scale = Eigen::VectorXd::Ones(dimension);
    Eigen::VectorXd column_scale = Eigen::VectorXd::Ones(dimension);
    Matrix scaled = matrix;
    for (int pass = 0; pass < 4; ++pass) {
        for (Eigen::Index row = 0; row < scaled.rows(); ++row) {
            const double norm = scaled.row(row).cwiseAbs().maxCoeff();
            if (norm > 0.0 && std::isfinite(norm)) {
                row_scale[row] /= norm;
                scaled.row(row) /= norm;
            }
        }
        for (Eigen::Index column = 0; column < scaled.cols(); ++column) {
            const double norm = scaled.col(column).cwiseAbs().maxCoeff();
            if (norm > 0.0 && std::isfinite(norm)) {
                column_scale[column] /= norm;
                scaled.col(column) /= norm;
            }
        }
    }
    Eigen::JacobiSVD<Matrix> decomposition(
        scaled, Eigen::ComputeFullU | Eigen::ComputeFullV
    );
    decomposition.setThreshold(
        4096.0 * std::numeric_limits<double>::epsilon()
        * static_cast<double>(dimension)
    );
    result.rank = static_cast<std::size_t>(decomposition.rank());
    const auto singular_values = decomposition.singularValues();
    if (result.rank == dimension && singular_values.size() != 0
        && singular_values[dimension - 1] > 0.0) {
        result.condition_number =
            singular_values[0] / singular_values[dimension - 1];
    }
    if (result.rank != dimension
        || !std::isfinite(result.condition_number)
        || result.condition_number > kSolvedStateSensitivityConditionNumberMax) {
        return result;
    }
    const Matrix scaled_rhs = (-row_scale).asDiagonal() * partials;
    const Matrix scaled_solution = decomposition.solve(scaled_rhs);
    const Matrix solution = column_scale.asDiagonal() * scaled_solution;
    const Matrix residual = matrix * solution + partials;
    const double numerator = residual.cwiseAbs().maxCoeff();
    const double denominator = std::max(
        1.0,
        matrix.cwiseAbs().maxCoeff() * solution.cwiseAbs().maxCoeff()
            + partials.cwiseAbs().maxCoeff()
    );
    result.backward_error = numerator / denominator;
    const double tolerance = solved_state_sensitivity_backward_error_max(
        dimension
    );
    if (!std::isfinite(result.backward_error)
        || result.backward_error > tolerance || !solution.allFinite()) {
        return result;
    }
    result.solution.assign(solution.data(), solution.data() + solution.size());
    return result;
}

}  // namespace

SolvedStateSensitivity evaluate_solved_state_sensitivity(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    const std::vector<std::string>& component_ids,
    const std::vector<double>& feed_amounts_mol,
    const std::vector<SolvedPhaseState>& phases
) {
    SolvedStateSensitivity result;
    const ProviderActiveParameterSet active = provider.active_parameters();
    result.parameter_identities = active.identities;
    if (active.identities.empty()) {
        result.failure_reason = "no_active_parameters";
        return result;
    }
    if (!std::isfinite(temperature_k) || temperature_k <= 0.0
        || !std::isfinite(pressure_pa) || pressure_pa <= 0.0
        || component_ids.size() != feed_amounts_mol.size()) {
        result.failure_reason = "invalid_solved_problem";
        return result;
    }
    if (!validate_phases(feed_amounts_mol, phases, result)) {
        return result;
    }
    const std::size_t component_count = component_ids.size();
    const std::size_t phase_count = phases.size();
    const std::size_t parameter_count = active.identities.size();
    if (phase_count > 1 && provider.transport().component_charges != nullptr
        && std::any_of(
            provider.transport().component_charges,
            provider.transport().component_charges + component_count,
            [](std::int32_t charge) { return charge != 0; }
        )) {
        result.failure_reason = "unsupported_charged_multiphase_constraints";
        return result;
    }
    result.reference_phase_index = phase_count - 1;
    result.dimension = (phase_count - 1) * component_count + phase_count;
    append_identities(result, component_ids, phase_count);
    std::vector<ParameterizedPhaseEvaluation> evaluations;
    evaluations.reserve(phase_count);
    try {
        for (const auto& phase : phases) {
            evaluations.push_back(provider.evaluate_parameterized_phase(
                temperature_k,
                phase.amounts_mol,
                phase.volume_m3,
                kPackingFractionMinimum,
                kPackingFractionMaximum,
                active
            ));
        }
    } catch (const std::exception&) {
        result.failure_reason = "provider_partial_evaluation_failed";
        return result;
    }
    if (!validate_defining_equations(
            evaluations, pressure_pa, component_count, result
        )) {
        return result;
    }
    result.residual_state_jacobian.assign(
        result.dimension * result.dimension, 0.0
    );
    result.residual_parameter_jacobian.assign(
        result.dimension * parameter_count, 0.0
    );
    const std::size_t reference = phase_count - 1;
    const std::size_t physical_count = component_count + 1;
    std::size_t residual = 0;
    for (std::size_t phase = 0; phase < reference; ++phase) {
        for (std::size_t component = 0; component < component_count; ++component) {
            for (std::size_t independent = 0; independent < reference; ++independent) {
                for (std::size_t varied = 0; varied < component_count; ++varied) {
                    double derivative = evaluations[reference].phase.hessian[
                        component * physical_count + varied
                    ];
                    if (independent == phase) {
                        derivative += evaluations[phase].phase.hessian[
                            component * physical_count + varied
                        ];
                    }
                    result.residual_state_jacobian[
                        residual * result.dimension
                        + amount_state_index(independent, varied, component_count)
                    ] = derivative;
                }
            }
            result.residual_state_jacobian[
                residual * result.dimension
                + volume_state_index(phase, phase_count, component_count)
            ] = evaluations[phase].phase.hessian[
                component * physical_count + component_count
            ];
            result.residual_state_jacobian[
                residual * result.dimension
                + volume_state_index(reference, phase_count, component_count)
            ] = -evaluations[reference].phase.hessian[
                component * physical_count + component_count
            ];
            for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
                result.residual_parameter_jacobian[
                    residual * parameter_count + parameter
                ] = evaluations[phase]
                        .chemical_potential_parameter_derivatives_over_rt[
                            component * parameter_count + parameter
                        ]
                    - evaluations[reference]
                        .chemical_potential_parameter_derivatives_over_rt[
                            component * parameter_count + parameter
                        ];
            }
            ++residual;
        }
    }
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        if (phase < reference) {
            for (std::size_t component = 0; component < component_count; ++component) {
                result.residual_state_jacobian[
                    residual * result.dimension
                    + amount_state_index(phase, component, component_count)
                ] = -kGasConstantJPerMolK * temperature_k
                    * evaluations[phase].phase.hessian[
                        component_count * physical_count + component
                    ] / pressure_pa;
            }
        } else {
            for (std::size_t independent = 0; independent < reference; ++independent) {
                for (std::size_t component = 0; component < component_count; ++component) {
                    result.residual_state_jacobian[
                        residual * result.dimension
                        + amount_state_index(independent, component, component_count)
                    ] = kGasConstantJPerMolK * temperature_k
                        * evaluations[reference].phase.hessian[
                            component_count * physical_count + component
                        ] / pressure_pa;
                }
            }
        }
        result.residual_state_jacobian[
            residual * result.dimension
            + volume_state_index(phase, phase_count, component_count)
        ] = -kGasConstantJPerMolK * temperature_k
            * evaluations[phase].phase.hessian[
                component_count * physical_count + component_count
            ] / pressure_pa;
        for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
            result.residual_parameter_jacobian[
                residual * parameter_count + parameter
            ] = evaluations[phase].pressure_parameter_derivatives_pa[parameter]
                / pressure_pa;
        }
        ++residual;
    }
    const RankRevealingSolve solved = solve_total_derivatives(
        result.residual_state_jacobian,
        result.residual_parameter_jacobian,
        result.dimension,
        parameter_count
    );
    result.rank = solved.rank;
    result.scaled_condition_number_2 = solved.condition_number;
    result.linear_solve_backward_error = solved.backward_error;
    if (solved.rank != result.dimension) {
        result.failure_reason = "singular_residual_jacobian";
        return result;
    }
    if (!std::isfinite(solved.condition_number)
        || solved.condition_number > kSolvedStateSensitivityConditionNumberMax) {
        result.failure_reason = "ill_conditioned_residual_jacobian";
        return result;
    }
    if (solved.solution.empty()) {
        result.failure_reason = "implicit_linear_solve_failed";
        return result;
    }
    result.state_parameter_jacobian = solved.solution;
    const std::size_t output_count = phase_count * (component_count + 1);
    result.output_parameter_jacobian.assign(output_count * parameter_count, 0.0);
    for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
        for (std::size_t phase = 0; phase < reference; ++phase) {
            for (std::size_t component = 0; component < component_count; ++component) {
                const double derivative = solved.solution[
                    amount_state_index(phase, component, component_count)
                        * parameter_count + parameter
                ];
                result.output_parameter_jacobian[
                    output_amount_index(phase, component, component_count)
                        * parameter_count + parameter
                ] = derivative;
                result.output_parameter_jacobian[
                    output_amount_index(reference, component, component_count)
                        * parameter_count + parameter
                ] -= derivative;
            }
        }
        for (std::size_t phase = 0; phase < phase_count; ++phase) {
            result.output_parameter_jacobian[
                output_volume_index(phase, component_count)
                    * parameter_count + parameter
            ] = solved.solution[
                volume_state_index(phase, phase_count, component_count)
                    * parameter_count + parameter
            ];
        }
    }
    result.status = "available";
    result.failure_reason.clear();
    return result;
}

}  // namespace epcsaft_engine::equilibrium
