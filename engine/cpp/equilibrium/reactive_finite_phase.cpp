#include "chemical_equilibrium.hpp"
#include "provider.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <optional>
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
    bool phase_stacked = false;
    std::vector<double> reference_potentials;
};

DenseMatrix admissible_invariant_matrix(
    const PhaseDirectionCompilation& directions,
    const std::vector<int>& charges
) {
    const std::size_t phase_count = directions.phase_ids.size();
    const std::size_t species_count = directions.species_ids.size();
    const std::size_t dimension = phase_count * species_count;
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
        invariants(
            directions.invariant_matrix.values.data(),
            static_cast<Eigen::Index>(directions.invariant_matrix.rows),
            static_cast<Eigen::Index>(dimension)
        );
    double charge_norm_squared = 0.0;
    for (int charge : charges) {
        charge_norm_squared += static_cast<double>(charge * charge);
    }
    if (charge_norm_squared == 0.0) {
        return directions.invariant_matrix;
    }
    Eigen::MatrixXd charge_basis = Eigen::MatrixXd::Zero(
        static_cast<Eigen::Index>(phase_count),
        static_cast<Eigen::Index>(dimension)
    );
    const double inverse_charge_norm = 1.0 / std::sqrt(charge_norm_squared);
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        for (std::size_t species = 0; species < species_count; ++species) {
            charge_basis(
                static_cast<Eigen::Index>(phase),
                static_cast<Eigen::Index>(phase * species_count + species)
            ) = static_cast<double>(charges[species]) * inverse_charge_norm;
        }
    }
    const Eigen::MatrixXd projected = invariants
        * (Eigen::MatrixXd::Identity(
            static_cast<Eigen::Index>(dimension),
            static_cast<Eigen::Index>(dimension)
        ) - charge_basis.transpose() * charge_basis);
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        projected, Eigen::ComputeFullV
    );
    const double threshold = 256.0 * std::numeric_limits<double>::epsilon()
        * static_cast<double>(dimension)
        * (svd.singularValues().size() == 0 ? 1.0 : svd.singularValues()(0));
    std::size_t rank = 0;
    for (Eigen::Index index = 0; index < svd.singularValues().size(); ++index) {
        rank += static_cast<std::size_t>(svd.singularValues()(index) > threshold);
    }
    Eigen::MatrixXd basis = svd.matrixV()
        .leftCols(static_cast<Eigen::Index>(rank))
        .transpose();
    DenseMatrix result;
    result.rows = rank;
    result.columns = dimension;
    result.values.reserve(rank * dimension);
    for (Eigen::Index row = 0; row < basis.rows(); ++row) {
        for (Eigen::Index column = 0; column < basis.cols(); ++column) {
            result.values.push_back(basis(row, column));
        }
    }
    return result;
}

NormalizedBalances normalized_balances(
    const CompiledReactionSystem& system,
    const PhaseDirectionCompilation* directions
) {
    NormalizedBalances result{
        directions == nullptr
            ? system.balance_matrix
            : admissible_invariant_matrix(*directions, system.charges),
        directions == nullptr
            ? system.balance_totals
            : std::vector<double>{},
        {},
        directions != nullptr,
        directions == nullptr ? system.g_ref : directions->reference_potentials,
    };
    if (directions != nullptr) {
        result.totals.assign(result.matrix.rows, 0.0);
        for (std::size_t row = 0; row < result.matrix.rows; ++row) {
            for (std::size_t coordinate = 0;
                 coordinate < result.matrix.columns;
                 ++coordinate) {
                result.totals[row] += result.matrix(row, coordinate)
                    * directions->phase_reference_amounts[coordinate];
            }
        }
    }
    result.scales.assign(result.matrix.rows, 1.0);
    for (std::size_t row = 0; row < result.matrix.rows; ++row) {
        double scale = std::abs(result.totals[row]);
        for (std::size_t column = 0; column < result.matrix.columns; ++column) {
            const double reference_amount = directions == nullptr
                ? system.feed_amounts[column]
                : directions->phase_reference_amounts[column];
            scale += std::abs(result.matrix(row, column) * reference_amount);
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

void shift_extensive_amount_coordinates(
    const AmountChart& chart,
    std::vector<double>& coordinates,
    std::size_t offset,
    double shift
) {
    if (!chart.ionic()) {
        for (std::size_t species = 0; species < chart.charges.size(); ++species) {
            coordinates[offset + species] += shift;
        }
        return;
    }
    coordinates[offset] += shift;
    const std::size_t neutral_offset = 1
        + chart.cation_indices.size() - 1
        + chart.anion_indices.size() - 1;
    for (std::size_t neutral = 0;
         neutral < chart.neutral_indices.size();
         ++neutral) {
        coordinates[offset + neutral_offset + neutral] += shift;
    }
}

struct PhaseBlock {
    std::vector<double> amounts;
    std::vector<double> amount_hessians;
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
    const AmountChart& chart,
    const std::vector<double>& g_ref,
    const std::vector<double>& g_ref_parameter_derivatives,
    double temperature_k,
    double pressure_pa,
    const std::array<double, 2>& packing_interval,
    const std::vector<double>& variables,
    std::size_t offset
) {
    const std::size_t species_count = g_ref.size();
    const std::size_t amount_dimension = chart.coordinate_count();
    const std::size_t phase_dimension = amount_dimension + 1;
    const std::size_t physical_dimension = species_count + 1;
    const std::size_t parameter_count = active_parameters.identities.size();
    if (!g_ref_parameter_derivatives.empty()
        && g_ref_parameter_derivatives.size() != species_count * parameter_count) {
        throw std::invalid_argument(
            "reactive finite-phase reference derivative dimensions changed"
        );
    }
    PhaseBlock result;
    const std::vector<double> amount_coordinates(
        variables.begin() + static_cast<std::ptrdiff_t>(offset),
        variables.begin() + static_cast<std::ptrdiff_t>(offset + amount_dimension)
    );
    const AmountChartEvaluation amount_chart =
        evaluate_amount_chart(chart, amount_coordinates);
    result.amounts = amount_chart.amounts;
    result.amount_hessians = amount_chart.amount_hessians;
    result.log_packing = variables[offset + amount_dimension];
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
    const std::size_t geometry_dimension = physical_dimension + parameter_count;
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
    if (mechanical.gradient.size() != physical_dimension
        || mechanical.hessian.size() != physical_dimension * physical_dimension
        || (parameter_count != 0
            && state_parameter_cross.size() != physical_dimension * parameter_count)) {
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

    std::vector<double> geometry_input_jacobian(
        physical_dimension * phase_dimension, 0.0
    );
    for (std::size_t species = 0; species < species_count; ++species) {
        for (std::size_t coordinate = 0;
             coordinate < amount_dimension;
             ++coordinate) {
            geometry_input_jacobian[species * phase_dimension + coordinate] =
                amount_chart.jacobian[species * amount_dimension + coordinate];
        }
    }
    geometry_input_jacobian[
        species_count * phase_dimension + amount_dimension
    ] = 1.0;

    std::vector<double> physical_from_chart(
        physical_dimension * phase_dimension, 0.0
    );
    for (std::size_t species = 0; species < species_count; ++species) {
        for (std::size_t coordinate = 0;
             coordinate < amount_dimension;
             ++coordinate) {
            physical_from_chart[species * phase_dimension + coordinate] =
                amount_chart.jacobian[species * amount_dimension + coordinate];
        }
    }
    for (std::size_t coordinate = 0; coordinate < phase_dimension; ++coordinate) {
        for (std::size_t physical = 0;
             physical < physical_dimension;
             ++physical) {
            physical_from_chart[species_count * phase_dimension + coordinate] +=
                geometry.gradient[physical]
                * geometry_input_jacobian[physical * phase_dimension + coordinate];
        }
    }

    std::vector<double> volume_hessian(
        phase_dimension * phase_dimension, 0.0
    );
    for (std::size_t left = 0; left < phase_dimension; ++left) {
        for (std::size_t right = 0; right < phase_dimension; ++right) {
            for (std::size_t physical_left = 0;
                 physical_left < physical_dimension;
                 ++physical_left) {
                for (std::size_t physical_right = 0;
                     physical_right < physical_dimension;
                     ++physical_right) {
                    volume_hessian[left * phase_dimension + right] +=
                        geometry_input_jacobian[
                            physical_left * phase_dimension + left
                        ]
                        * geometry.hessian[
                            physical_left * geometry_dimension + physical_right
                        ]
                        * geometry_input_jacobian[
                            physical_right * phase_dimension + right
                        ];
                }
            }
            if (left < amount_dimension && right < amount_dimension) {
                for (std::size_t species = 0;
                     species < species_count;
                     ++species) {
                    volume_hessian[left * phase_dimension + right] +=
                        geometry.gradient[species]
                        * amount_chart.amount_hessians[
                            species * amount_dimension * amount_dimension
                                + left * amount_dimension + right
                        ];
                }
            }
        }
    }
    result.gradient.resize(phase_dimension, 0.0);
    result.hessian.resize(phase_dimension * phase_dimension, 0.0);
    result.state_jacobian = physical_from_chart;
    for (std::size_t coordinate = 0; coordinate < phase_dimension; ++coordinate) {
        for (std::size_t physical = 0;
             physical < physical_dimension;
             ++physical) {
            result.gradient[coordinate] +=
                physical_from_chart[physical * phase_dimension + coordinate]
                * result.physical_gradient[physical];
        }
        for (std::size_t right = 0; right < phase_dimension; ++right) {
            for (std::size_t physical_left = 0;
                 physical_left < physical_dimension;
                 ++physical_left) {
                for (std::size_t physical_right = 0;
                     physical_right < physical_dimension;
                     ++physical_right) {
                    result.hessian[coordinate * phase_dimension + right] +=
                        physical_from_chart[
                            physical_left * phase_dimension + coordinate
                        ]
                        * mechanical.hessian[
                            physical_left * physical_dimension + physical_right
                        ]
                        * physical_from_chart[
                            physical_right * phase_dimension + right
                        ];
                }
            }
            if (coordinate < amount_dimension && right < amount_dimension) {
                for (std::size_t species = 0;
                     species < species_count;
                     ++species) {
                    result.hessian[coordinate * phase_dimension + right] +=
                        result.physical_gradient[species]
                        * amount_chart.amount_hessians[
                            species * amount_dimension * amount_dimension
                                + coordinate * amount_dimension + right
                        ];
                }
            }
            result.hessian[coordinate * phase_dimension + right] +=
                result.physical_gradient.back()
                * volume_hessian[coordinate * phase_dimension + right];
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
        const std::size_t parameter_coordinate = physical_dimension + parameter;
        const double volume_direct = geometry.gradient[parameter_coordinate];
        result.volume_parameter_direct[parameter] = volume_direct;
        std::vector<double> physical_gradient_parameter(
            physical_dimension, 0.0
        );
        for (std::size_t physical = 0;
             physical < physical_dimension;
             ++physical) {
            physical_gradient_parameter[physical] =
                mechanical.hessian[
                    physical * physical_dimension + species_count
                ]
                    * volume_direct
                + state_parameter_cross[physical * parameter_count + parameter];
        }
        for (std::size_t coordinate = 0;
             coordinate < phase_dimension;
             ++coordinate) {
            double volume_jacobian_parameter_derivative = 0.0;
            for (std::size_t physical = 0;
                 physical < physical_dimension;
                 ++physical) {
                volume_jacobian_parameter_derivative +=
                    geometry.hessian[
                        physical * geometry_dimension + parameter_coordinate
                    ]
                    * geometry_input_jacobian[
                        physical * phase_dimension + coordinate
                    ];
            }
            double derivative = volume_jacobian_parameter_derivative
                * result.physical_gradient.back();
            for (std::size_t physical = 0;
                 physical < physical_dimension;
                 ++physical) {
                derivative +=
                    physical_from_chart[physical * phase_dimension + coordinate]
                    * physical_gradient_parameter[physical];
            }
            if (!g_ref_parameter_derivatives.empty()) {
                for (std::size_t species = 0;
                     species < species_count;
                     ++species) {
                    derivative += physical_from_chart[
                        species * phase_dimension + coordinate
                    ] * g_ref_parameter_derivatives[
                        species * parameter_count + parameter
                    ];
                }
            }
            result.objective_parameter_cross[
                coordinate * parameter_count + parameter
            ] = derivative;
        }
        if (!g_ref_parameter_derivatives.empty()) {
            for (std::size_t species = 0;
                 species < species_count;
                 ++species) {
                result.chemical_potential_parameter_direct[
                    species * parameter_count + parameter
                ] += g_ref_parameter_derivatives[
                    species * parameter_count + parameter
                ];
            }
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
    const AmountChart chart = make_amount_chart(system.charges);
    const std::size_t phase_dimension = chart.coordinate_count() + 1;
    const std::size_t variable_count = phase_count * phase_dimension;
    const std::size_t constraint_count = balances.matrix.rows;
    if (variables.size() != variable_count
        || multipliers.size() != constraint_count) {
        throw std::invalid_argument("reactive finite-phase NLP dimensions changed");
    }
    MultiphaseEvaluation result;
    result.gradient.assign(variable_count, 0.0);
    result.constraints.assign(constraint_count, 0.0);
    result.jacobian.assign(constraint_count * variable_count, 0.0);
    result.objective_hessian.assign(variable_count * variable_count, 0.0);
    result.phases.reserve(phase_count);
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        const std::size_t offset = phase * phase_dimension;
        std::vector<double> phase_reference = system.g_ref;
        std::vector<double> phase_reference_derivatives =
            system.g_ref_parameter_derivatives;
        if (balances.phase_stacked) {
            phase_reference.assign(system.species_count, 0.0);
            phase_reference_derivatives.clear();
            for (std::size_t species = 0; species < system.species_count; ++species) {
                phase_reference[species] = balances.reference_potentials[
                    phase * system.species_count + species
                ];
            }
        }
        result.phases.push_back(evaluate_phase(
            provider,
            active_parameters,
            chart,
            phase_reference,
            phase_reference_derivatives,
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
                const std::size_t balance_column = balances.phase_stacked
                    ? phase * species_count + species
                    : species;
                result.constraints[row] += balances.matrix(row, balance_column)
                    * block.amounts[species];
                for (std::size_t coordinate = 0;
                     coordinate + 1 < phase_dimension;
                     ++coordinate) {
                    result.jacobian[
                        row * variable_count + offset + coordinate
                    ] += balances.matrix(row, balance_column)
                        * block.state_jacobian[
                            species * phase_dimension + coordinate
                        ];
                }
            }
        }
    }
    result.lagrangian_hessian = result.objective_hessian;
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        const std::size_t offset = phase * phase_dimension;
        const std::size_t amount_dimension = phase_dimension - 1;
        for (std::size_t left = 0; left < amount_dimension; ++left) {
            for (std::size_t right = 0; right < amount_dimension; ++right) {
                for (std::size_t row = 0; row < balances.matrix.rows; ++row) {
                    for (std::size_t species = 0;
                         species < species_count;
                         ++species) {
                        result.lagrangian_hessian[
                            (offset + left) * variable_count + offset + right
                        ] += multipliers[row]
                            * balances.matrix(
                                row,
                                balances.phase_stacked
                                    ? phase * species_count + species
                                    : species
                            )
                            * result.phases[phase].amount_hessians[
                                species * amount_dimension * amount_dimension
                                    + left * amount_dimension + right
                            ];
                    }
                }
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
        for (std::size_t row = 0; row < baseline.constraints.size(); ++row) {
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
            const MultiphaseEvaluation& evaluated = evaluate(x, zero_multipliers());
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
            const MultiphaseEvaluation& evaluated = evaluate(x, zero_multipliers());
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
            const MultiphaseEvaluation& evaluated = evaluate(x, zero_multipliers());
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
            const MultiphaseEvaluation& evaluated = evaluate(x, lambda);
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

    const MultiphaseEvaluation& evaluate(
        const Ipopt::Number* x,
        const std::vector<double>& multipliers
    ) const {
        const std::vector<double> point(x, x + initial_.size());
        const bool zero = std::all_of(
            multipliers.begin(), multipliers.end(), [](double value) {
                return value == 0.0;
            }
        );
        std::vector<double>& cached_point = zero ? base_point_ : lagrangian_point_;
        std::vector<double>& cached_multipliers =
            zero ? base_multipliers_ : lagrangian_multipliers_;
        std::optional<MultiphaseEvaluation>& cached =
            zero ? base_evaluation_ : lagrangian_evaluation_;
        if (!cached.has_value() || cached_point != point
            || cached_multipliers != multipliers) {
            cached = evaluate_multiphase(
                provider_,
                active_parameters_,
                system_,
                balances_,
                temperature_k_,
                pressure_pa_,
                packing_intervals_,
                point,
                multipliers,
                phase_count_
            );
            cached_point = point;
            cached_multipliers = multipliers;
        }
        return *cached;
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
    mutable std::vector<double> base_point_;
    mutable std::vector<double> base_multipliers_;
    mutable std::optional<MultiphaseEvaluation> base_evaluation_;
    mutable std::vector<double> lagrangian_point_;
    mutable std::vector<double> lagrangian_multipliers_;
    mutable std::optional<MultiphaseEvaluation> lagrangian_evaluation_;
};

double phase_distinctness(const std::vector<PhaseBlock>& phases) {
    if (phases.size() < 2) {
        return std::numeric_limits<double>::infinity();
    }
    double result = std::numeric_limits<double>::infinity();
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
            result = std::min(result, distance);
        }
    }
    return result;
}

std::vector<double> equality_multipliers_at(
    const MultiphaseEvaluation& evaluation
) {
    const std::size_t constraint_count = evaluation.constraints.size();
    const std::size_t variable_count = evaluation.gradient.size();
    if (constraint_count == 0) {
        return {};
    }
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
        jacobian(
            evaluation.jacobian.data(),
            static_cast<Eigen::Index>(constraint_count),
            static_cast<Eigen::Index>(variable_count)
        );
    Eigen::Map<const Eigen::VectorXd> gradient(
        evaluation.gradient.data(), static_cast<Eigen::Index>(variable_count)
    );
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> factor(
        jacobian.transpose()
    );
    if (static_cast<std::size_t>(factor.rank()) != constraint_count) {
        throw std::runtime_error(
            "reactive finite-phase equality multiplier recovery is rank deficient"
        );
    }
    const Eigen::VectorXd recovered = factor.solve(-gradient);
    if (!recovered.allFinite()) {
        throw std::runtime_error(
            "reactive finite-phase equality multiplier recovery is nonfinite"
        );
    }
    return std::vector<double>(recovered.data(), recovered.data() + recovered.size());
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
    const std::size_t physical_dimension = species_count + 1;
    const std::size_t phase_dimension =
        evaluation.gradient.size() / phase_count;
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
            std::vector<double> physical_response(physical_dimension, 0.0);
            for (std::size_t physical = 0;
                 physical < physical_dimension;
                 ++physical) {
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
                for (std::size_t physical = 0;
                     physical < physical_dimension;
                     ++physical) {
                    chemical_response +=
                        block.physical_hessian[
                            species * physical_dimension + physical
                        ]
                        * physical_response[physical];
                }
                result.chemical_potential_derivatives[
                    (phase * species_count + species) * parameter_count + parameter
                ] = chemical_response;
            }
            result.volume_derivatives[phase * parameter_count + parameter] =
                physical_response.back();
            double pressure_response = block.pressure_parameter_direct[parameter];
            for (std::size_t physical = 0;
                 physical < physical_dimension;
                 ++physical) {
                pressure_response -= kGasConstantJPerMolK * temperature_k
                    * block.physical_hessian[
                        species_count * physical_dimension + physical
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
    bool perform_derivative_audits,
    const PhaseDirectionCompilation* phase_directions
) {
    ReactiveFinitePhaseResult result;
    const std::size_t phase_count = starts.size();
    const std::size_t species_count = system.species_count;
    const AmountChart chart = make_amount_chart(system.charges);
    const std::size_t amount_dimension = chart.coordinate_count();
    const std::size_t phase_dimension = amount_dimension + 1;
    const std::size_t variable_count = phase_count * phase_dimension;
    if (phase_count < 2 || packing_fraction_intervals.size() != phase_count) {
        throw std::invalid_argument(
            "reactive finite-phase solve requires at least two declared phases"
        );
    }
    if (system.species_count != system.original_species_count
        || !system.removed_species_indices.empty()) {
        throw std::invalid_argument(
            "reactive finite-phase solve requires a complete strict-interior species face"
        );
    }
    if (phase_directions != nullptr
        && (phase_directions->phase_ids.size() != phase_count
            || phase_directions->species_ids.size() != species_count
            || phase_directions->charges.size() != species_count
            || phase_directions->phase_reference_amounts.size()
                != phase_count * species_count
            || phase_directions->reference_potentials.size()
                != phase_count * species_count)) {
        throw std::invalid_argument(
            "reactive finite-phase direction dimensions changed"
        );
    }
    if (!std::isfinite(temperature_k) || temperature_k <= 0.0
        || !std::isfinite(pressure_pa) || pressure_pa <= 0.0
        || !std::isfinite(trace_floor) || trace_floor <= 0.0
        || max_iterations <= 0) {
        throw std::invalid_argument("reactive finite-phase state controls are invalid");
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
    const double total_feed_amount_mol = std::accumulate(
        system.feed_amounts.begin(), system.feed_amounts.end(), 0.0
    );
    if (!std::isfinite(total_feed_mass_kg) || total_feed_mass_kg <= 0.0
        || !std::isfinite(total_feed_amount_mol) || total_feed_amount_mol <= 0.0) {
        throw std::invalid_argument("reactive finite-phase feed amount is invalid");
    }
    const double amount_log_scale = std::log(total_feed_amount_mol);
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
        const std::vector<double> amount_coordinates =
            invert_amount_chart(chart, start.amounts_mol);
        std::copy(
            amount_coordinates.begin(),
            amount_coordinates.end(),
            initial.begin() + static_cast<std::ptrdiff_t>(offset)
        );
        std::fill(
            lower.begin() + static_cast<std::ptrdiff_t>(offset),
            lower.begin() + static_cast<std::ptrdiff_t>(offset + amount_dimension),
            -40.0
        );
        std::fill(
            upper.begin() + static_cast<std::ptrdiff_t>(offset),
            upper.begin() + static_cast<std::ptrdiff_t>(offset + amount_dimension),
            40.0
        );
        if (!chart.ionic()) {
            for (std::size_t species = 0; species < species_count; ++species) {
                lower[offset + species] = std::log(0.1 * trace_floor);
                upper[offset + species] = std::log(std::nextafter(
                    total_feed_mass_kg / system.molar_masses_kg_per_mol[species],
                    std::numeric_limits<double>::infinity()
                ));
            }
        } else {
            double charge_equivalent_upper = 0.0;
            for (std::size_t species : chart.cation_indices) {
                charge_equivalent_upper +=
                    static_cast<double>(system.charges[species])
                    * total_feed_mass_kg
                    / system.molar_masses_kg_per_mol[species];
            }
            if (!std::isfinite(charge_equivalent_upper)
                || charge_equivalent_upper <= 0.0) {
                throw std::invalid_argument(
                    "ionic charge-equivalent upper bound is invalid"
                );
            }
            lower[offset] = std::log(0.1 * trace_floor);
            upper[offset] = std::log(charge_equivalent_upper);
            const std::size_t neutral_offset =
                1 + chart.cation_indices.size() - 1
                + chart.anion_indices.size() - 1;
            for (std::size_t neutral = 0;
                 neutral < chart.neutral_indices.size();
                 ++neutral) {
                const std::size_t species = chart.neutral_indices[neutral];
                lower[offset + neutral_offset + neutral] =
                    std::log(0.1 * trace_floor);
                upper[offset + neutral_offset + neutral] =
                    std::log(std::nextafter(
                        total_feed_mass_kg
                            / system.molar_masses_kg_per_mol[species],
                        std::numeric_limits<double>::infinity()
                    ));
            }
        }
        shift_extensive_amount_coordinates(
            chart, initial, offset, -amount_log_scale
        );
        shift_extensive_amount_coordinates(
            chart, lower, offset, -amount_log_scale
        );
        shift_extensive_amount_coordinates(
            chart, upper, offset, -amount_log_scale
        );
        initial[offset + amount_dimension] = std::log(start_packing);
        lower[offset + amount_dimension] = std::log(packing[0]);
        upper[offset + amount_dimension] = std::nextafter(
            std::log(packing[1]), -std::numeric_limits<double>::infinity()
        );
    }

    NormalizedBalances balances = normalized_balances(system, phase_directions);
    for (double& coefficient : balances.matrix.values) {
        coefficient *= total_feed_amount_mol;
    }
    const std::size_t constraint_count = balances.matrix.rows;
    if (phase_count * chart.coordinate_count() <= constraint_count) {
        throw std::invalid_argument(
            "configured reactive phases have no independent composition directions"
        );
    }
    result.admissible_nullity = phase_count * chart.coordinate_count()
        - constraint_count;
    result.declared_direction_rank = phase_directions == nullptr
        ? result.admissible_nullity
        : phase_directions->direction_rank;
    const MultiphaseEvaluation initial_evaluation = evaluate_multiphase(
        provider,
        parameters,
        system,
        balances,
        temperature_k,
        pressure_pa,
        packing_fraction_intervals,
        initial,
        std::vector<double>(constraint_count, 0.0),
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
    application->Options()->SetStringValue("nlp_scaling_method", "none");
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

    std::vector<double> solution = raw_problem->solution();
    std::vector<double> multipliers = raw_problem->multipliers();
    // Fixed-pressure phase amounts do not depend on the packing coordinate.
    // Polish each phase's scalar packing coordinate after IPOPT so the direct
    // pressure check is as accurate as the otherwise well-scaled KKT solve.
    for (int iteration = 0; iteration < 8; ++iteration) {
        const MultiphaseEvaluation pressure_evaluation = evaluate_multiphase(
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
        double maximum_relative_residual = 0.0;
        bool changed = false;
        for (std::size_t phase = 0; phase < phase_count; ++phase) {
            const PhaseBlock& block = pressure_evaluation.phases[phase];
            const double residual = block.pressure - pressure_pa;
            maximum_relative_residual = std::max(
                maximum_relative_residual, std::abs(residual) / pressure_pa
            );
            const double volume_derivative = block.state_jacobian[
                species_count * phase_dimension + amount_dimension
            ];
            const double pressure_derivative =
                -kGasConstantJPerMolK * temperature_k
                * block.physical_hessian[
                    species_count * (species_count + 1) + species_count
                ]
                * volume_derivative;
            if (!std::isfinite(pressure_derivative)
                || std::abs(pressure_derivative)
                    <= std::numeric_limits<double>::min()) {
                continue;
            }
            const std::size_t coordinate = phase * phase_dimension + amount_dimension;
            const double unconstrained_step = -residual / pressure_derivative;
            const double lower_margin = solution[coordinate] - lower[coordinate];
            const double upper_margin = upper[coordinate] - solution[coordinate];
            const double step = std::clamp(
                unconstrained_step,
                -0.5 * lower_margin,
                0.5 * upper_margin
            );
            if (std::isfinite(step) && step != 0.0) {
                solution[coordinate] += step;
                changed = true;
            }
        }
        if (maximum_relative_residual <= 0.1 * kChemicalPressureTolerance
            || !changed) {
            break;
        }
    }
    const MultiphaseEvaluation polished_primal = evaluate_multiphase(
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
    multipliers = equality_multipliers_at(polished_primal);
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
    result.objective_over_rt = total_feed_amount_mol * evaluated.objective;
    result.balance_inf_norm = inf_norm(evaluated.constraints);
    result.phase_distinctness = phase_distinctness(evaluated.phases);
    result.phase_amounts_mol.reserve(phase_count);
    result.phase_volumes_m3.reserve(phase_count);
    result.phase_pressures_pa.reserve(phase_count);
    result.phase_packing_fractions.reserve(phase_count);
    result.phase_chemical_potentials_over_rt.reserve(phase_count);
    std::vector<double> stationarity = evaluated.gradient;
    for (std::size_t row = 0; row < evaluated.constraints.size(); ++row) {
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
        std::vector<double> physical_amounts = block.amounts;
        for (double& amount : physical_amounts) {
            amount *= total_feed_amount_mol;
        }
        result.phase_amounts_mol.push_back(std::move(physical_amounts));
        result.phase_volumes_m3.push_back(
            total_feed_amount_mol * block.volume
        );
        result.phase_pressures_pa.push_back(block.pressure);
        result.phase_packing_fractions.push_back(block.packing);
        result.phase_chemical_potentials_over_rt.push_back(block.chemical_potentials);
        result.pressure_relative_inf_norm = std::max(
            result.pressure_relative_inf_norm,
            std::abs(block.pressure - pressure_pa) / pressure_pa
        );
        for (double amount : block.amounts) {
            minimum_amount = std::min(
                minimum_amount, total_feed_amount_mol * amount
            );
        }
        const std::size_t offset = phase * phase_dimension;
        for (std::size_t coordinate = 0; coordinate < phase_dimension; ++coordinate) {
            minimum_bound_margin = std::min({
                minimum_bound_margin,
                solution[offset + coordinate] - lower[offset + coordinate],
                upper[offset + coordinate] - solution[offset + coordinate],
            });
        }
        if (phase_directions == nullptr) {
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
    }
    if (phase_directions != nullptr) {
        for (std::size_t direction = 0;
             direction < phase_directions->direction_matrix.rows;
             ++direction) {
            double affinity = 0.0;
            for (std::size_t phase = 0; phase < phase_count; ++phase) {
                for (std::size_t species = 0; species < species_count; ++species) {
                    affinity += phase_directions->direction_matrix(
                        direction, phase * species_count + species
                    ) * evaluated.phases[phase].chemical_potentials[species];
                }
            }
            result.reaction_affinity_inf_norm = std::max(
                result.reaction_affinity_inf_norm, std::abs(affinity)
            );
        }
        result.transfer_chemical_potential_inf_norm =
            result.reaction_affinity_inf_norm;
    }
    const double charge_norm_squared = std::inner_product(
        system.charges.begin(), system.charges.end(), system.charges.begin(), 0.0
    );
    for (std::size_t phase = 1;
         phase_directions == nullptr && phase < phase_count;
         ++phase) {
        double electric_offset = 0.0;
        if (charge_norm_squared > 0.0) {
            for (std::size_t species = 0; species < species_count; ++species) {
                electric_offset += static_cast<double>(system.charges[species])
                    * (evaluated.phases[phase].chemical_potentials[species]
                        - evaluated.phases[0].chemical_potentials[species]);
            }
            electric_offset /= charge_norm_squared;
        }
        for (std::size_t species = 0; species < species_count; ++species) {
            const double residual =
                evaluated.phases[phase].chemical_potentials[species]
                - evaluated.phases[0].chemical_potentials[species]
                - static_cast<double>(system.charges[species]) * electric_offset;
            result.transfer_chemical_potential_inf_norm = std::max(
                result.transfer_chemical_potential_inf_norm,
                std::abs(residual)
            );
        }
    }
    result.reduced_hessian_min_eigenvalue = reduced_hessian_minimum(
        evaluated, variable_count, evaluated.constraints.size()
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
    for (double& derivative : result.amount_derivatives) {
        derivative *= total_feed_amount_mol;
    }
    for (double& derivative : result.volume_derivatives) {
        derivative *= total_feed_amount_mol;
    }
    if (!parameters.identities.empty() && result.sensitivity_status != "available") {
        result.accepted = false;
        result.failure_reason = "reactive_finite_phase_sensitivity_failed";
    }
    result.provider_evaluation_count = static_cast<std::size_t>(provider.evaluation_count());
    return result;
}

}  // namespace epcsaft_engine::equilibrium
