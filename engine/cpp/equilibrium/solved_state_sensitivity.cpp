#include "solved_state_sensitivity.hpp"

#include "held2.hpp"
#include "phase_equilibrium.hpp"
#include "physical_constants.hpp"

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

using epcsaft::detail::GAS_CONSTANT;
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
    const std::vector<epcsaft::NativeParameterizedPhaseState>& evaluations,
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

std::vector<double> held2_independent(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical
) {
    const auto modified = held2_transform_physical_fractions(
        coordinates, physical
    );
    std::vector<double> result;
    result.reserve(coordinates.independent_indices.size());
    for (std::size_t eos : coordinates.independent_indices) {
        const auto retained = std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            eos
        );
        result.push_back(modified[static_cast<std::size_t>(
            retained - coordinates.retained_indices.begin()
        )]);
    }
    return result;
}

std::vector<double> modified_potential_rows(
    const Held2Coordinates& coordinates
) {
    const std::size_t component_count = coordinates.charges.size();
    std::vector<double> result(
        coordinates.retained_indices.size() * component_count, 0.0
    );
    const double eliminated_charge =
        coordinates.charges[coordinates.eliminated_index];
    for (std::size_t row = 0;
         row < coordinates.retained_indices.size(); ++row) {
        const std::size_t retained = coordinates.retained_indices[row];
        result[row * component_count + retained] =
            1.0 / coordinates.modified_factors[row];
        result[row * component_count + coordinates.eliminated_index] =
            -coordinates.charges[retained]
            / eliminated_charge
            / coordinates.modified_factors[row];
    }
    return result;
}

SolvedStateSensitivity evaluate_charged_sensitivity(
    const epcsaft::NativeModel& eos,
    double temperature_k,
    double pressure_pa,
    const std::vector<std::string>& component_ids,
    const std::vector<double>& feed_amounts_mol,
    const std::vector<SolvedPhaseState>& phases,
    const std::vector<epcsaft::NativeActiveCoordinate>& active
) {
    SolvedStateSensitivity result;
    result.parameter_identities = epcsaft::active_parameter_identities(active);
    const std::size_t component_count = component_ids.size();
    const std::size_t phase_count = phases.size();
    const std::size_t parameter_count = active.size();
    const std::vector<double> charges = phase_eos_component_charges(eos);
    Held2Coordinates coordinates;
    try {
        coordinates = make_held2_coordinates(charges, component_ids);
    } catch (const std::exception&) {
        result.failure_reason = "unsupported_charged_coordinate_chart";
        return result;
    }
    const std::size_t independent_count =
        coordinates.independent_indices.size();
    const std::size_t potential_count =
        coordinates.retained_indices.size();
    const std::size_t block_size = independent_count + 2;
    result.dimension = phase_count * block_size;
    result.reference_phase_index = phase_count - 1;
    const double feed_total = std::accumulate(
        feed_amounts_mol.begin(), feed_amounts_mol.end(), 0.0
    );
    const auto composition_jacobian =
        held2_physical_fraction_jacobian(coordinates);
    const auto potential_rows = modified_potential_rows(coordinates);
    std::vector<double> phase_amounts;
    std::vector<double> molar_volumes;
    std::vector<std::vector<double>> compositions;
    std::vector<std::vector<double>> independent_values;
    std::vector<epcsaft::NativeParameterizedPhaseState> evaluations;
    try {
        for (const auto& phase : phases) {
            const double amount = std::accumulate(
                phase.amounts_mol.begin(), phase.amounts_mol.end(), 0.0
            );
            std::vector<double> composition = phase.amounts_mol;
            for (double& value : composition) {
                value /= amount;
            }
            const auto chart = held2_independent(coordinates, composition);
            static_cast<void>(held2_lift_independent_fractions(
                coordinates, chart
            ));
            const double molar_volume = phase.volume_m3 / amount;
            phase_amounts.push_back(amount);
            molar_volumes.push_back(molar_volume);
            compositions.push_back(composition);
            independent_values.push_back(chart);
            evaluations.push_back(eos.evaluate_parameterized_phase_state(
                temperature_k,
                composition,
                molar_volume,
                active
            ));
            result.eos_evaluation_count += evaluations.back().evaluations.total();
        }
    } catch (const std::exception&) {
        result.failure_reason = "eos_partial_evaluation_failed";
        return result;
    }
    std::vector<std::vector<double>> modified_potentials;
    for (const auto& evaluation : evaluations) {
        modified_potentials.push_back(held2_transform_modified_potentials(
            coordinates,
            std::vector<double>(
                evaluation.phase.gradient.begin(),
                evaluation.phase.gradient.begin()
                    + static_cast<std::ptrdiff_t>(component_count)
            )
        ));
    }
    for (std::size_t phase = 0; phase + 1 < phase_count; ++phase) {
        for (std::size_t potential = 0;
             potential < potential_count; ++potential) {
            if (std::abs(
                    modified_potentials[phase][potential]
                    - modified_potentials.back()[potential]
                ) > kChemicalPotentialResidualMax) {
                result.failure_reason = "defining_equations_not_satisfied";
                return result;
            }
        }
    }
    for (const auto& evaluation : evaluations) {
        if (std::abs(evaluation.phase.pressure_pa / pressure_pa - 1.0)
            > kPressureResidualMax) {
            result.failure_reason = "defining_equations_not_satisfied";
            return result;
        }
    }

    result.residual_state_jacobian.assign(
        result.dimension * result.dimension, 0.0
    );
    result.residual_parameter_jacobian.assign(
        result.dimension * parameter_count, 0.0
    );
    auto beta_index = [block_size](std::size_t phase) {
        return phase * block_size;
    };
    auto composition_index = [block_size](
        std::size_t phase, std::size_t coordinate
    ) {
        return phase * block_size + 1 + coordinate;
    };
    auto log_volume_index = [block_size, independent_count](
        std::size_t phase
    ) {
        return phase * block_size + 1 + independent_count;
    };
    std::size_t residual = 0;
    result.residual_identities.push_back("total_phase_amount_balance");
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        result.residual_state_jacobian[
            residual * result.dimension + beta_index(phase)
        ] = 1.0;
    }
    ++residual;
    const auto feed_composition = [&] {
        std::vector<double> value = feed_amounts_mol;
        for (double& item : value) {
            item /= feed_total;
        }
        return value;
    }();
    const auto feed_independent = held2_independent(
        coordinates, feed_composition
    );
    if (std::abs(
            std::accumulate(
                phase_amounts.begin(), phase_amounts.end(), 0.0
            ) - feed_total
        ) > 2.0e-10 * std::max(1.0, std::abs(feed_total))) {
        result.failure_reason = "material_balance_mismatch";
        return result;
    }
    for (std::size_t coordinate = 0;
         coordinate < independent_count; ++coordinate) {
        result.residual_identities.push_back(
            "modified_material_balance["
            + component_ids[coordinates.independent_indices[coordinate]]
            + "]"
        );
        for (std::size_t phase = 0; phase < phase_count; ++phase) {
            result.residual_state_jacobian[
                residual * result.dimension + beta_index(phase)
            ] = independent_values[phase][coordinate];
            result.residual_state_jacobian[
                residual * result.dimension
                + composition_index(phase, coordinate)
            ] = phase_amounts[phase];
        }
        double balance = -feed_total * feed_independent[coordinate];
        for (std::size_t phase = 0; phase < phase_count; ++phase) {
            balance += phase_amounts[phase]
                * independent_values[phase][coordinate];
        }
        if (std::abs(balance)
            > 2.0e-10 * std::max(1.0, std::abs(feed_total))) {
            result.failure_reason = "material_balance_mismatch";
            return result;
        }
        ++residual;
    }
    const std::size_t reference = phase_count - 1;
    for (std::size_t phase = 0; phase < reference; ++phase) {
        for (std::size_t potential = 0;
             potential < potential_count; ++potential) {
            result.residual_identities.push_back(
                "phase[" + std::to_string(phase)
                + "].modified_mu["
                + component_ids[coordinates.retained_indices[potential]]
                + "]-phase[" + std::to_string(reference)
                + "].modified_mu"
            );
            for (const std::size_t affected : {phase, reference}) {
                const double sign = affected == phase ? 1.0 : -1.0;
                const std::size_t physical_count = component_count + 1;
                for (std::size_t coordinate = 0;
                     coordinate < independent_count; ++coordinate) {
                    double derivative = 0.0;
                    for (std::size_t raw = 0; raw < component_count; ++raw) {
                        const double transform = potential_rows[
                            potential * component_count + raw
                        ];
                        for (std::size_t component = 0;
                             component < component_count; ++component) {
                            derivative += transform
                                * evaluations[affected].phase.hessian[
                                    raw * physical_count + component
                                ]
                                * composition_jacobian[
                                    component * independent_count + coordinate
                                ];
                        }
                    }
                    result.residual_state_jacobian[
                        residual * result.dimension
                        + composition_index(affected, coordinate)
                    ] += sign * derivative;
                }
                double volume_derivative = 0.0;
                for (std::size_t raw = 0; raw < component_count; ++raw) {
                    volume_derivative += potential_rows[
                        potential * component_count + raw
                    ] * evaluations[affected].phase.hessian[
                        raw * physical_count + component_count
                    ] * molar_volumes[affected];
                }
                result.residual_state_jacobian[
                    residual * result.dimension
                    + log_volume_index(affected)
                ] += sign * volume_derivative;
                for (std::size_t parameter = 0;
                     parameter < parameter_count; ++parameter) {
                    double derivative = 0.0;
                    for (std::size_t raw = 0;
                         raw < component_count; ++raw) {
                        derivative += potential_rows[
                            potential * component_count + raw
                        ] * evaluations[affected]
                            .chemical_potential_parameter_derivatives_over_rt[
                                raw * parameter_count + parameter
                            ];
                    }
                    result.residual_parameter_jacobian[
                        residual * parameter_count + parameter
                    ] += sign * derivative;
                }
            }
            ++residual;
        }
    }
    const double rt = GAS_CONSTANT * temperature_k;
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        result.residual_identities.push_back(
            "phase[" + std::to_string(phase)
            + "].pressure_over_target_minus_one"
        );
        const std::size_t physical_count = component_count + 1;
        for (std::size_t coordinate = 0;
             coordinate < independent_count; ++coordinate) {
            double derivative = 0.0;
            for (std::size_t component = 0;
                 component < component_count; ++component) {
                derivative += -rt * evaluations[phase].phase.hessian[
                    component_count * physical_count + component
                ] * composition_jacobian[
                    component * independent_count + coordinate
                ] / pressure_pa;
            }
            result.residual_state_jacobian[
                residual * result.dimension
                + composition_index(phase, coordinate)
            ] = derivative;
        }
        result.residual_state_jacobian[
            residual * result.dimension + log_volume_index(phase)
        ] = -rt * evaluations[phase].phase.hessian[
            component_count * physical_count + component_count
        ] * molar_volumes[phase] / pressure_pa;
        for (std::size_t parameter = 0;
             parameter < parameter_count; ++parameter) {
            result.residual_parameter_jacobian[
                residual * parameter_count + parameter
            ] = evaluations[phase]
                .pressure_parameter_derivatives_pa[parameter] / pressure_pa;
        }
        ++residual;
    }
    if (residual != result.dimension) {
        result.failure_reason = "charged_residual_layout_incomplete";
        return result;
    }
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        result.state_identities.push_back(
            "phase[" + std::to_string(phase) + "].amount_mol"
        );
        for (std::size_t coordinate = 0;
             coordinate < independent_count; ++coordinate) {
            result.state_identities.push_back(
                "phase[" + std::to_string(phase) + "].modified_fraction["
                + component_ids[coordinates.independent_indices[coordinate]]
                + "]"
            );
        }
        result.state_identities.push_back(
            "phase[" + std::to_string(phase) + "].log_molar_volume"
        );
        for (const auto& component : component_ids) {
            result.output_identities.push_back(
                "phase[" + std::to_string(phase) + "].amount["
                + component + "]"
            );
        }
        result.output_identities.push_back(
            "phase[" + std::to_string(phase) + "].volume_m3"
        );
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
    result.output_parameter_jacobian.assign(
        output_count * parameter_count, 0.0
    );
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        for (std::size_t parameter = 0;
             parameter < parameter_count; ++parameter) {
            for (std::size_t component = 0;
                 component < component_count; ++component) {
                double derivative = compositions[phase][component]
                    * solved.solution[
                        beta_index(phase) * parameter_count + parameter
                    ];
                for (std::size_t coordinate = 0;
                     coordinate < independent_count; ++coordinate) {
                    derivative += phase_amounts[phase]
                        * composition_jacobian[
                            component * independent_count + coordinate
                        ] * solved.solution[
                            composition_index(phase, coordinate)
                                * parameter_count + parameter
                        ];
                }
                result.output_parameter_jacobian[
                    output_amount_index(
                        phase, component, component_count
                    ) * parameter_count + parameter
                ] = derivative;
            }
            result.output_parameter_jacobian[
                output_volume_index(phase, component_count)
                    * parameter_count + parameter
            ] = molar_volumes[phase] * solved.solution[
                    beta_index(phase) * parameter_count + parameter
                ] + phase_amounts[phase] * molar_volumes[phase]
                    * solved.solution[
                        log_volume_index(phase) * parameter_count + parameter
                    ];
        }
    }
    result.status = "available";
    result.failure_reason.clear();
    return result;
}

}  // namespace

bool solve_sensitivity_linear_system(
    SolvedStateSensitivity& result,
    std::size_t parameter_count
) {
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
        return false;
    }
    if (!std::isfinite(solved.condition_number)
        || solved.condition_number > kSolvedStateSensitivityConditionNumberMax) {
        result.failure_reason = "ill_conditioned_residual_jacobian";
        return false;
    }
    if (solved.solution.empty()) {
        result.failure_reason = "implicit_linear_solve_failed";
        return false;
    }
    result.state_parameter_jacobian = solved.solution;
    return true;
}

static SolvedStateSensitivity evaluate_state_sensitivity(
    const epcsaft::NativeModel& eos,
    double temperature_k,
    double pressure_pa,
    const std::vector<std::string>& component_ids,
    const std::vector<double>& feed_amounts_mol,
    const std::vector<SolvedPhaseState>& phases,
    bool allow_charged_multiphase
) {
    SolvedStateSensitivity result;
    const std::vector<epcsaft::NativeActiveCoordinate> active = eos.active_coordinates();
    result.parameter_identities = epcsaft::active_parameter_identities(active);
    if (active.empty()) {
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
    const std::size_t parameter_count = active.size();
    if (phase_count > 1
        && std::any_of(
            eos.charge_numbers().begin(),
            eos.charge_numbers().end(),
            [](int charge) { return charge != 0; }
        )) {
        if (!allow_charged_multiphase) {
            result.failure_reason =
                "unsupported_charged_multiphase_constraints";
            return result;
        }
        return evaluate_charged_sensitivity(
            eos,
            temperature_k,
            pressure_pa,
            component_ids,
            feed_amounts_mol,
            phases,
            active
        );
    }
    result.reference_phase_index = phase_count - 1;
    result.dimension = (phase_count - 1) * component_count + phase_count;
    append_identities(result, component_ids, phase_count);
    std::vector<epcsaft::NativeParameterizedPhaseState> evaluations;
    evaluations.reserve(phase_count);
    try {
        for (const auto& phase : phases) {
            evaluations.push_back(eos.evaluate_parameterized_phase_state(
                temperature_k,
                phase.amounts_mol,
                phase.volume_m3,
                active
            ));
            result.eos_evaluation_count += evaluations.back().evaluations.total();
        }
    } catch (const std::exception&) {
        result.failure_reason = "eos_partial_evaluation_failed";
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
                ] = -GAS_CONSTANT * temperature_k
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
                    ] = GAS_CONSTANT * temperature_k
                        * evaluations[reference].phase.hessian[
                            component_count * physical_count + component
                        ] / pressure_pa;
                }
            }
        }
        result.residual_state_jacobian[
            residual * result.dimension
            + volume_state_index(phase, phase_count, component_count)
        ] = -GAS_CONSTANT * temperature_k
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

SolvedStateSensitivity evaluate_solved_state_sensitivity(
    const epcsaft::NativeModel& eos,
    double temperature_k,
    double pressure_pa,
    const std::vector<std::string>& component_ids,
    const std::vector<double>& feed_amounts_mol,
    const std::vector<SolvedPhaseState>& phases
) {
    return evaluate_state_sensitivity(
        eos,
        temperature_k,
        pressure_pa,
        component_ids,
        feed_amounts_mol,
        phases,
        false
    );
}

SolvedStateSensitivity evaluate_declared_state_sensitivity(
    const epcsaft::NativeModel& eos,
    double temperature_k,
    double pressure_pa,
    const std::vector<std::string>& component_ids,
    const std::vector<double>& feed_amounts_mol,
    const std::vector<SolvedPhaseState>& phases
) {
    return evaluate_state_sensitivity(
        eos,
        temperature_k,
        pressure_pa,
        component_ids,
        feed_amounts_mol,
        phases,
        true
    );
}

}  // namespace epcsaft_engine::equilibrium
