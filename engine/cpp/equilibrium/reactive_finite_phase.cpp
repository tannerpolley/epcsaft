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
#include <Highs.h>
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpTNLP.hpp>

namespace epcsaft_engine::equilibrium {
namespace {

constexpr double kGasConstantJPerMolK = 8.31446261815324;
constexpr double kReactiveFinitePhaseStartBalanceTolerance = 1.0e-8;
constexpr double kReactiveFinitePhaseSolverTolerance = 1.0e-10;
constexpr double kReactiveFinitePhaseAcceptableTolerance = 1.0e-8;
constexpr double kReactiveFinitePhaseDistinctnessMinimum = 1.0e-6;
constexpr double kReactiveFinitePhaseBoundMarginMinimum = 1.0e-7;
constexpr double kReactiveFinitePhaseStrictMinimumFactor = 1.0e-10;

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

std::optional<std::vector<double>> strict_positive_seed(
    const CompiledReactionSystem& system,
    const std::vector<std::size_t>& support,
    double trace_floor,
    double ion_fraction_max
) {
    const std::size_t supported_count = support.size();
    const std::size_t minimum_column = supported_count;
    const bool ion_limited = std::isfinite(ion_fraction_max);
    const std::size_t row_count = system.balance_matrix.rows
        + supported_count + (ion_limited ? 1 : 0);
    HighsModel model;
    model.lp_.num_col_ = static_cast<HighsInt>(supported_count + 1);
    model.lp_.num_row_ = static_cast<HighsInt>(row_count);
    model.lp_.sense_ = ObjSense::kMaximize;
    model.lp_.col_cost_.assign(supported_count + 1, 0.0);
    model.lp_.col_cost_[minimum_column] = 1.0;
    model.lp_.col_lower_.assign(supported_count + 1, 0.0);
    model.lp_.col_upper_.assign(supported_count + 1, kHighsInf);
    model.lp_.row_lower_ = system.balance_totals;
    model.lp_.row_upper_ = system.balance_totals;
    model.lp_.row_lower_.insert(
        model.lp_.row_lower_.end(), supported_count, 0.0
    );
    model.lp_.row_upper_.insert(
        model.lp_.row_upper_.end(), supported_count, kHighsInf
    );
    if (ion_limited) {
        model.lp_.row_lower_.push_back(-kHighsInf);
        model.lp_.row_upper_.push_back(0.0);
    }
    model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
    model.lp_.a_matrix_.start_ = {0};
    for (std::size_t column = 0; column < supported_count; ++column) {
        const std::size_t species = support[column];
        for (std::size_t row = 0; row < system.balance_matrix.rows; ++row) {
            const double value = system.balance_matrix(row, species);
            if (value != 0.0) {
                model.lp_.a_matrix_.index_.push_back(static_cast<HighsInt>(row));
                model.lp_.a_matrix_.value_.push_back(value);
            }
        }
        model.lp_.a_matrix_.index_.push_back(static_cast<HighsInt>(
            system.balance_matrix.rows + column
        ));
        model.lp_.a_matrix_.value_.push_back(1.0);
        if (ion_limited) {
            model.lp_.a_matrix_.index_.push_back(static_cast<HighsInt>(row_count - 1));
            model.lp_.a_matrix_.value_.push_back(
                system.charges[species] == 0
                    ? -ion_fraction_max
                    : 1.0 - ion_fraction_max
            );
        }
        model.lp_.a_matrix_.start_.push_back(static_cast<HighsInt>(
            model.lp_.a_matrix_.index_.size()
        ));
    }
    for (std::size_t column = 0; column < supported_count; ++column) {
        model.lp_.a_matrix_.index_.push_back(static_cast<HighsInt>(
            system.balance_matrix.rows + column
        ));
        model.lp_.a_matrix_.value_.push_back(-1.0);
    }
    model.lp_.a_matrix_.start_.push_back(static_cast<HighsInt>(
        model.lp_.a_matrix_.index_.size()
    ));
    Highs highs;
    if (highs.setOptionValue("output_flag", false) == HighsStatus::kError
        || highs.setOptionValue("threads", 1) == HighsStatus::kError
        || highs.setOptionValue("solver", std::string("simplex"))
            == HighsStatus::kError
        || highs.passModel(model) == HighsStatus::kError
        || highs.run() == HighsStatus::kError
        || highs.getModelStatus() != HighsModelStatus::kOptimal) {
        return std::nullopt;
    }
    const std::vector<double>& solution = highs.getSolution().col_value;
    if (solution.size() != supported_count + 1
        || !std::isfinite(solution.back()) || solution.back() <= trace_floor) {
        return std::nullopt;
    }
    std::vector<double> amounts(system.species_count, 0.0);
    for (std::size_t column = 0; column < supported_count; ++column) {
        if (!std::isfinite(solution[column]) || solution[column] <= trace_floor) {
            return std::nullopt;
        }
        amounts[support[column]] = solution[column];
    }
    return amounts;
}

DenseMatrix admissible_invariant_matrix(
    const PhaseDirectionCompilation& directions,
    const std::vector<int>& charges,
    const std::vector<double>& molar_masses_kg_per_mol,
    const std::vector<std::vector<std::size_t>>& supports
) {
    const std::size_t phase_count = directions.phase_ids.size();
    const std::size_t species_count = directions.species_ids.size();
    const std::size_t dimension = phase_count * species_count;
    std::vector<std::size_t> active_coordinates;
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        for (std::size_t species : supports[phase]) {
            active_coordinates.push_back(phase * species_count + species);
        }
    }
    Eigen::MatrixXd invariants(
        static_cast<Eigen::Index>(directions.invariant_matrix.rows),
        static_cast<Eigen::Index>(active_coordinates.size())
    );
    for (std::size_t row = 0; row < directions.invariant_matrix.rows; ++row) {
        for (std::size_t column = 0; column < active_coordinates.size(); ++column) {
            invariants(
                static_cast<Eigen::Index>(row),
                static_cast<Eigen::Index>(column)
            ) = directions.invariant_matrix(row, active_coordinates[column]);
        }
    }
    Eigen::MatrixXd charge_basis = Eigen::MatrixXd::Zero(
        static_cast<Eigen::Index>(phase_count),
        static_cast<Eigen::Index>(active_coordinates.size())
    );
    std::size_t active_offset = 0;
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        double charge_norm_squared = 0.0;
        for (std::size_t species : supports[phase]) {
            charge_norm_squared += static_cast<double>(charges[species] * charges[species]);
        }
        const double inverse_charge_norm = charge_norm_squared == 0.0
            ? 0.0
            : 1.0 / std::sqrt(charge_norm_squared);
        for (std::size_t local = 0; local < supports[phase].size(); ++local) {
            charge_basis(
                static_cast<Eigen::Index>(phase),
                static_cast<Eigen::Index>(active_offset + local)
            ) = static_cast<double>(charges[supports[phase][local]])
                * inverse_charge_norm;
        }
        active_offset += supports[phase].size();
    }
    const Eigen::MatrixXd charge_projector = Eigen::MatrixXd::Identity(
            static_cast<Eigen::Index>(active_coordinates.size()),
            static_cast<Eigen::Index>(active_coordinates.size())
        ) - charge_basis.transpose() * charge_basis;
    const Eigen::MatrixXd projected = invariants * charge_projector;
    Eigen::RowVectorXd mass(active_coordinates.size());
    for (std::size_t column = 0; column < active_coordinates.size(); ++column) {
        mass(static_cast<Eigen::Index>(column)) = molar_masses_kg_per_mol[
            active_coordinates[column] % species_count
        ];
    }
    const Eigen::RowVectorXd projected_mass = mass * charge_projector;
    if (projected_mass.squaredNorm() == 0.0) {
        throw std::invalid_argument(
            "configured reactive phases have no admissible positive mass invariant"
        );
    }
    const Eigen::MatrixXd mass_complement = Eigen::MatrixXd::Identity(
        static_cast<Eigen::Index>(active_coordinates.size()),
        static_cast<Eigen::Index>(active_coordinates.size())
    ) - projected_mass.transpose() * projected_mass
        / projected_mass.squaredNorm();
    const Eigen::MatrixXd remaining = projected * mass_complement;
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        remaining, Eigen::ComputeFullV
    );
    const double threshold = 256.0 * std::numeric_limits<double>::epsilon()
        * static_cast<double>(active_coordinates.size())
        * (svd.singularValues().size() == 0 ? 1.0 : svd.singularValues()(0));
    std::size_t rank = 0;
    for (Eigen::Index index = 0; index < svd.singularValues().size(); ++index) {
        rank += static_cast<std::size_t>(svd.singularValues()(index) > threshold);
    }
    Eigen::MatrixXd basis(static_cast<Eigen::Index>(rank + 1), mass.size());
    basis.row(0) = mass;
    if (rank != 0) {
        basis.bottomRows(static_cast<Eigen::Index>(rank)) = svd.matrixV()
            .leftCols(static_cast<Eigen::Index>(rank))
            .transpose();
    }
    DenseMatrix result;
    result.rows = rank + 1;
    result.columns = dimension;
    result.values.assign((rank + 1) * dimension, 0.0);
    for (Eigen::Index row = 0; row < basis.rows(); ++row) {
        for (Eigen::Index column = 0; column < basis.cols(); ++column) {
            result(
                static_cast<std::size_t>(row),
                active_coordinates[static_cast<std::size_t>(column)]
            ) = basis(row, column);
        }
    }
    return result;
}

NormalizedBalances normalized_balances(
    const CompiledReactionSystem& system,
    const PhaseDirectionCompilation* directions,
    const std::vector<std::vector<std::size_t>>& supports
) {
    NormalizedBalances result{
        directions == nullptr
            ? system.balance_matrix
            : admissible_invariant_matrix(
                *directions,
                system.charges,
                system.molar_masses_kg_per_mol,
                supports
            ),
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

struct PhaseLayout {
    std::vector<std::size_t> supported_species;
    AmountChart chart;
    std::size_t offset = 0;

    [[nodiscard]] std::size_t amount_dimension() const {
        return chart.coordinate_count();
    }

    [[nodiscard]] std::size_t dimension() const {
        return amount_dimension() + 1;
    }
};

std::vector<PhaseLayout> make_phase_layouts(
    const std::vector<int>& charges,
    const std::vector<std::vector<std::size_t>>& supports
) {
    std::vector<PhaseLayout> layouts;
    layouts.reserve(supports.size());
    std::size_t offset = 0;
    for (const std::vector<std::size_t>& support : supports) {
        if (support.empty()) {
            throw std::invalid_argument("reactive finite-phase support must not be empty");
        }
        std::vector<int> supported_charges;
        supported_charges.reserve(support.size());
        std::size_t previous = charges.size();
        for (std::size_t species : support) {
            if (species >= charges.size() || (previous != charges.size() && species <= previous)) {
                throw std::invalid_argument(
                    "reactive finite-phase support indices must be unique Provider-order indices"
                );
            }
            previous = species;
            supported_charges.push_back(charges[species]);
        }
        PhaseLayout layout{support, make_amount_chart(supported_charges), offset};
        offset += layout.dimension();
        layouts.push_back(std::move(layout));
    }
    return layouts;
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
    std::vector<double> objective_pressure_cross;
    std::vector<double> objective_parameter_cross;
    std::vector<double> volume_parameter_direct;
    std::vector<double> pressure_parameter_direct;
    std::vector<double> chemical_potential_parameter_direct;
};

PhaseBlock evaluate_phase(
    const ProviderContext& provider,
    const ProviderActiveParameterSet& active_parameters,
    const PhaseLayout& layout,
    const std::vector<double>& g_ref,
    const std::vector<double>& g_ref_pressure_derivatives_per_pa,
    const std::vector<double>& g_ref_parameter_derivatives,
    double temperature_k,
    double pressure_pa,
    const std::array<double, 2>& packing_interval,
    const std::vector<double>& variables,
    std::size_t offset
) {
    const std::size_t species_count = g_ref.size();
    const AmountChart& chart = layout.chart;
    const std::vector<std::size_t>& supported_species = layout.supported_species;
    const std::size_t support_count = supported_species.size();
    const std::size_t amount_dimension = chart.coordinate_count();
    const std::size_t phase_dimension = amount_dimension + 1;
    const std::size_t physical_dimension = support_count + 1;
    const std::size_t parameter_count = active_parameters.identities.size();
    if (!g_ref_parameter_derivatives.empty()
        && g_ref_parameter_derivatives.size() != species_count * parameter_count) {
        throw std::invalid_argument(
            "reactive finite-phase reference derivative dimensions changed"
        );
    }
    if (!g_ref_pressure_derivatives_per_pa.empty()
        && g_ref_pressure_derivatives_per_pa.size() != species_count) {
        throw std::invalid_argument(
            "reactive finite-phase reference pressure derivative dimensions changed"
        );
    }
    PhaseBlock result;
    const std::vector<double> amount_coordinates(
        variables.begin() + static_cast<std::ptrdiff_t>(offset),
        variables.begin() + static_cast<std::ptrdiff_t>(offset + amount_dimension)
    );
    const AmountChartEvaluation amount_chart =
        evaluate_amount_chart(chart, amount_coordinates);
    result.amounts.assign(species_count, 0.0);
    result.amount_hessians.assign(
        species_count * amount_dimension * amount_dimension, 0.0
    );
    for (std::size_t supported = 0; supported < support_count; ++supported) {
        const std::size_t species = supported_species[supported];
        result.amounts[species] = amount_chart.amounts[supported];
        std::copy_n(
            amount_chart.amount_hessians.begin()
                + static_cast<std::ptrdiff_t>(
                    supported * amount_dimension * amount_dimension
                ),
            amount_dimension * amount_dimension,
            result.amount_hessians.begin()
                + static_cast<std::ptrdiff_t>(
                    species * amount_dimension * amount_dimension
                )
        );
    }
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
        || geometry.hessian.size() != geometry_dimension * geometry_dimension
        || geometry.supported_component_indices != supported_species) {
        throw std::invalid_argument(
            "reactive finite-phase inverse-packing tensor dimensions changed"
        );
    }
    result.volume = geometry.volume_m3;

    const ParameterizedPhaseEvaluation parameterized =
        provider.evaluate_parameterized_phase(
            temperature_k,
            result.amounts,
            result.volume,
            packing_interval[0],
            packing_interval[1],
            active_parameters
        );
    if (parameterized.supported_component_indices != supported_species) {
        throw std::invalid_argument("reactive finite-phase Provider support changed");
    }
    MixturePhaseEvaluation mechanical = parameterized.phase;
    const std::vector<double>& state_parameter_cross =
        parameterized.state_parameter_derivatives;
    result.pressure_parameter_direct =
        parameterized.pressure_parameter_derivatives_pa;
    result.chemical_potential_parameter_direct.assign(
        species_count * parameter_count, 0.0
    );
    for (std::size_t supported = 0; supported < support_count; ++supported) {
        for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
            result.chemical_potential_parameter_direct[
                supported_species[supported] * parameter_count + parameter
            ] = parameterized.chemical_potential_parameter_derivatives_over_rt[
                supported * parameter_count + parameter
            ];
        }
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
    for (std::size_t supported = 0; supported < support_count; ++supported) {
        result.physical_gradient[supported] += g_ref[supported_species[supported]];
    }
    result.physical_gradient.back() += pressure_over_rt;
    result.objective = mechanical.value + pressure_over_rt * result.volume;
    for (std::size_t supported = 0; supported < support_count; ++supported) {
        const std::size_t species = supported_species[supported];
        result.objective += g_ref[species] * result.amounts[species];
    }

    std::vector<double> geometry_input_jacobian(
        physical_dimension * phase_dimension, 0.0
    );
    for (std::size_t supported = 0; supported < support_count; ++supported) {
        for (std::size_t coordinate = 0;
             coordinate < amount_dimension;
             ++coordinate) {
            geometry_input_jacobian[supported * phase_dimension + coordinate] =
                amount_chart.jacobian[supported * amount_dimension + coordinate];
        }
    }
    geometry_input_jacobian[
        support_count * phase_dimension + amount_dimension
    ] = 1.0;

    std::vector<double> physical_from_chart(
        physical_dimension * phase_dimension, 0.0
    );
    for (std::size_t supported = 0; supported < support_count; ++supported) {
        for (std::size_t coordinate = 0;
             coordinate < amount_dimension;
             ++coordinate) {
            physical_from_chart[supported * phase_dimension + coordinate] =
                amount_chart.jacobian[supported * amount_dimension + coordinate];
        }
    }
    for (std::size_t coordinate = 0; coordinate < phase_dimension; ++coordinate) {
        for (std::size_t physical = 0;
             physical < physical_dimension;
             ++physical) {
            physical_from_chart[support_count * phase_dimension + coordinate] +=
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
                for (std::size_t supported = 0;
                     supported < support_count;
                     ++supported) {
                    volume_hessian[left * phase_dimension + right] +=
                        geometry.gradient[supported]
                        * amount_chart.amount_hessians[
                            supported * amount_dimension * amount_dimension
                                + left * amount_dimension + right
                        ];
                }
            }
        }
    }
    result.gradient.resize(phase_dimension, 0.0);
    result.hessian.resize(phase_dimension * phase_dimension, 0.0);
    result.state_jacobian.assign(
        (species_count + 1) * phase_dimension, 0.0
    );
    for (std::size_t supported = 0; supported < support_count; ++supported) {
        const std::size_t species = supported_species[supported];
        std::copy_n(
            physical_from_chart.begin()
                + static_cast<std::ptrdiff_t>(supported * phase_dimension),
            phase_dimension,
            result.state_jacobian.begin()
                + static_cast<std::ptrdiff_t>(species * phase_dimension)
        );
    }
    std::copy_n(
        physical_from_chart.begin()
            + static_cast<std::ptrdiff_t>(support_count * phase_dimension),
        phase_dimension,
        result.state_jacobian.begin()
            + static_cast<std::ptrdiff_t>(species_count * phase_dimension)
    );
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
                for (std::size_t supported = 0;
                     supported < support_count;
                     ++supported) {
                    result.hessian[coordinate * phase_dimension + right] +=
                        result.physical_gradient[supported]
                        * amount_chart.amount_hessians[
                            supported * amount_dimension * amount_dimension
                                + coordinate * amount_dimension + right
                        ];
                }
            }
            result.hessian[coordinate * phase_dimension + right] +=
                result.physical_gradient.back()
                * volume_hessian[coordinate * phase_dimension + right];
        }
    }

    result.chemical_potentials.assign(
        species_count, std::numeric_limits<double>::quiet_NaN()
    );
    for (std::size_t supported = 0; supported < support_count; ++supported) {
        const std::size_t species = supported_species[supported];
        result.chemical_potentials[species] =
            mechanical.gradient[supported] + g_ref[species];
    }

    result.objective_parameter_cross.assign(
        phase_dimension * parameter_count, 0.0
    );
    result.objective_pressure_cross.assign(phase_dimension, 0.0);
    for (std::size_t coordinate = 0; coordinate < phase_dimension; ++coordinate) {
        result.objective_pressure_cross[coordinate] =
            result.state_jacobian[species_count * phase_dimension + coordinate]
            / (kGasConstantJPerMolK * temperature_k);
        if (!g_ref_pressure_derivatives_per_pa.empty()) {
            for (std::size_t species = 0; species < species_count; ++species) {
                result.objective_pressure_cross[coordinate] +=
                    result.state_jacobian[species * phase_dimension + coordinate]
                    * g_ref_pressure_derivatives_per_pa[species];
            }
        }
    }
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
                    physical * physical_dimension + support_count
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
                for (std::size_t supported = 0;
                     supported < support_count;
                     ++supported) {
                    const std::size_t species = supported_species[supported];
                    derivative += physical_from_chart[
                        supported * phase_dimension + coordinate
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
            for (std::size_t supported = 0;
                 supported < support_count;
                 ++supported) {
                const std::size_t species = supported_species[supported];
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

PhaseBlock evaluate_ideal_phase(
    const PhaseLayout& layout,
    const std::vector<double>& g_ref,
    const std::vector<double>& g_ref_pressure_derivatives_per_pa,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& variables,
    std::size_t offset
) {
    const std::size_t species_count = g_ref.size();
    const std::size_t support_count = layout.supported_species.size();
    const std::size_t amount_dimension = layout.amount_dimension();
    const std::size_t phase_dimension = layout.dimension();
    const std::size_t physical_dimension = support_count + 1;
    const std::vector<double> amount_coordinates(
        variables.begin() + static_cast<std::ptrdiff_t>(offset),
        variables.begin() + static_cast<std::ptrdiff_t>(offset + amount_dimension)
    );
    const AmountChartEvaluation chart = evaluate_amount_chart(
        layout.chart, amount_coordinates
    );
    PhaseBlock result;
    result.amounts.assign(species_count, 0.0);
    result.amount_hessians.assign(
        species_count * amount_dimension * amount_dimension, 0.0
    );
    for (std::size_t supported = 0; supported < support_count; ++supported) {
        const std::size_t species = layout.supported_species[supported];
        result.amounts[species] = chart.amounts[supported];
        std::copy_n(
            chart.amount_hessians.begin()
                + static_cast<std::ptrdiff_t>(
                    supported * amount_dimension * amount_dimension
                ),
            amount_dimension * amount_dimension,
            result.amount_hessians.begin()
                + static_cast<std::ptrdiff_t>(
                    species * amount_dimension * amount_dimension
                )
        );
    }
    result.log_packing = variables[offset + amount_dimension];
    result.volume = std::exp(result.log_packing);
    if (!std::isfinite(result.volume) || result.volume <= 0.0) {
        throw std::domain_error("reactive ideal-gas volume chart is invalid");
    }
    result.packing = 0.0;
    const double total = std::accumulate(
        result.amounts.begin(), result.amounts.end(), 0.0
    );
    const double pressure_over_rt = pressure_pa
        / (kGasConstantJPerMolK * temperature_k);
    result.pressure = kGasConstantJPerMolK * temperature_k
        * total / result.volume;
    result.physical_gradient.assign(physical_dimension, 0.0);
    result.physical_hessian.assign(
        physical_dimension * physical_dimension, 0.0
    );
    result.chemical_potentials.assign(
        species_count, std::numeric_limits<double>::quiet_NaN()
    );
    result.objective = pressure_over_rt * result.volume;
    for (std::size_t supported = 0; supported < support_count; ++supported) {
        const std::size_t species = layout.supported_species[supported];
        const double amount = result.amounts[species];
        const double logarithm = std::log(amount / result.volume);
        result.objective += amount * (logarithm - 1.0 + g_ref[species]);
        result.physical_gradient[supported] = logarithm + g_ref[species];
        result.chemical_potentials[species] = logarithm + g_ref[species];
        result.physical_hessian[
            supported * physical_dimension + supported
        ] = 1.0 / amount;
        result.physical_hessian[
            supported * physical_dimension + support_count
        ] = -1.0 / result.volume;
        result.physical_hessian[
            support_count * physical_dimension + supported
        ] = -1.0 / result.volume;
    }
    result.physical_gradient.back() = -total / result.volume + pressure_over_rt;
    result.physical_hessian.back() = total / (result.volume * result.volume);

    result.state_jacobian.assign(
        (species_count + 1) * phase_dimension, 0.0
    );
    for (std::size_t supported = 0; supported < support_count; ++supported) {
        const std::size_t species = layout.supported_species[supported];
        std::copy_n(
            chart.jacobian.begin()
                + static_cast<std::ptrdiff_t>(supported * amount_dimension),
            amount_dimension,
            result.state_jacobian.begin()
                + static_cast<std::ptrdiff_t>(species * phase_dimension)
        );
    }
    result.state_jacobian[
        species_count * phase_dimension + amount_dimension
    ] = result.volume;
    result.gradient.assign(phase_dimension, 0.0);
    result.hessian.assign(phase_dimension * phase_dimension, 0.0);
    for (std::size_t left = 0; left < phase_dimension; ++left) {
        for (std::size_t physical = 0; physical < physical_dimension; ++physical) {
            const std::size_t full_physical = physical == support_count
                ? species_count
                : layout.supported_species[physical];
            result.gradient[left] += result.state_jacobian[
                full_physical * phase_dimension + left
            ] * result.physical_gradient[physical];
        }
        for (std::size_t right = 0; right < phase_dimension; ++right) {
            for (std::size_t physical_left = 0;
                 physical_left < physical_dimension;
                 ++physical_left) {
                const std::size_t full_left = physical_left == support_count
                    ? species_count
                    : layout.supported_species[physical_left];
                for (std::size_t physical_right = 0;
                     physical_right < physical_dimension;
                     ++physical_right) {
                    const std::size_t full_right = physical_right == support_count
                        ? species_count
                        : layout.supported_species[physical_right];
                    result.hessian[left * phase_dimension + right] +=
                        result.state_jacobian[
                            full_left * phase_dimension + left
                        ] * result.physical_hessian[
                            physical_left * physical_dimension + physical_right
                        ] * result.state_jacobian[
                            full_right * phase_dimension + right
                        ];
                }
            }
            if (left < amount_dimension && right < amount_dimension) {
                for (std::size_t supported = 0;
                     supported < support_count;
                     ++supported) {
                    result.hessian[left * phase_dimension + right] +=
                        result.physical_gradient[supported]
                        * chart.amount_hessians[
                            supported * amount_dimension * amount_dimension
                                + left * amount_dimension + right
                        ];
                }
            }
            if (left == amount_dimension && right == amount_dimension) {
                result.hessian[left * phase_dimension + right] +=
                    result.physical_gradient.back() * result.volume;
            }
        }
    }
    result.objective_pressure_cross.assign(phase_dimension, 0.0);
    for (std::size_t coordinate = 0; coordinate < phase_dimension; ++coordinate) {
        result.objective_pressure_cross[coordinate] = result.state_jacobian[
            species_count * phase_dimension + coordinate
        ] / (kGasConstantJPerMolK * temperature_k);
        if (!g_ref_pressure_derivatives_per_pa.empty()) {
            for (std::size_t species = 0; species < species_count; ++species) {
                result.objective_pressure_cross[coordinate] +=
                    result.state_jacobian[species * phase_dimension + coordinate]
                    * g_ref_pressure_derivatives_per_pa[species];
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
    const std::vector<PhaseLayout>& layouts
) {
    const std::size_t species_count = system.species_count;
    const std::size_t phase_count = layouts.size();
    const std::size_t variable_count = layouts.empty()
        ? 0
        : layouts.back().offset + layouts.back().dimension();
    const std::size_t constraint_count = balances.matrix.rows;
    if (variables.size() != variable_count
        || multipliers.size() != constraint_count
        || packing_intervals.size() != phase_count) {
        throw std::invalid_argument("reactive finite-phase NLP dimensions changed");
    }
    MultiphaseEvaluation result;
    result.gradient.assign(variable_count, 0.0);
    result.constraints.assign(constraint_count, 0.0);
    result.jacobian.assign(constraint_count * variable_count, 0.0);
    result.objective_hessian.assign(variable_count * variable_count, 0.0);
    result.phases.reserve(phase_count);
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        const PhaseLayout& layout = layouts[phase];
        const std::size_t offset = layout.offset;
        const std::size_t phase_dimension = layout.dimension();
        std::vector<double> phase_reference = system.g_ref;
        std::vector<double> phase_reference_pressure_derivatives =
            system.g_ref_pressure_derivatives_per_pa;
        std::vector<double> phase_reference_derivatives =
            system.g_ref_parameter_derivatives;
        if (balances.phase_stacked) {
            phase_reference.assign(system.species_count, 0.0);
            phase_reference_pressure_derivatives.clear();
            phase_reference_derivatives.clear();
            for (std::size_t species = 0; species < system.species_count; ++species) {
                phase_reference[species] = balances.reference_potentials[
                    phase * system.species_count + species
                ];
            }
        }
        result.phases.push_back(system.ideal_gas
            ? evaluate_ideal_phase(
                layout,
                phase_reference,
                phase_reference_pressure_derivatives,
                temperature_k,
                pressure_pa,
                variables,
                offset
            )
            : evaluate_phase(
                provider,
                active_parameters,
                layout,
                phase_reference,
                phase_reference_pressure_derivatives,
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
            const PhaseLayout& layout = layouts[phase];
            const std::size_t offset = layout.offset;
            const std::size_t phase_dimension = layout.dimension();
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
        const PhaseLayout& layout = layouts[phase];
        const std::size_t offset = layout.offset;
        const std::size_t amount_dimension = layout.amount_dimension();
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

void polish_single_phase_reaction_root(
    const ProviderContext& provider,
    const ProviderActiveParameterSet& active_parameters,
    const CompiledReactionSystem& system,
    const NormalizedBalances& balances,
    double temperature_k,
    double pressure_pa,
    const std::vector<std::array<double, 2>>& packing_intervals,
    const std::vector<PhaseLayout>& layouts,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    std::vector<double>& variables
) {
    if (layouts.size() != 1
        || balances.phase_stacked
        || balances.matrix.rows + system.reaction_matrix.rows + 1
            != variables.size()) {
        return;
    }
    const PhaseLayout& layout = layouts.front();
    const std::size_t dimension = variables.size();
    const std::size_t balance_count = balances.matrix.rows;
    const std::size_t reaction_count = system.reaction_matrix.rows;
    auto residual_and_jacobian = [&](const std::vector<double>& point) {
        const MultiphaseEvaluation evaluated = evaluate_multiphase(
            provider,
            active_parameters,
            system,
            balances,
            temperature_k,
            pressure_pa,
            packing_intervals,
            point,
            std::vector<double>(balance_count, 0.0),
            layouts
        );
        const PhaseBlock& phase = evaluated.phases.front();
        Eigen::VectorXd residual = Eigen::VectorXd::Zero(
            static_cast<Eigen::Index>(dimension)
        );
        Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(
            static_cast<Eigen::Index>(dimension),
            static_cast<Eigen::Index>(dimension)
        );
        for (std::size_t row = 0; row < balance_count; ++row) {
            residual(static_cast<Eigen::Index>(row)) = evaluated.constraints[row];
            for (std::size_t column = 0; column < dimension; ++column) {
                jacobian(
                    static_cast<Eigen::Index>(row),
                    static_cast<Eigen::Index>(column)
                ) = evaluated.jacobian[row * dimension + column];
            }
        }
        const std::size_t support_count = layout.supported_species.size();
        const std::size_t physical_dimension = support_count + 1;
        for (std::size_t reaction = 0; reaction < reaction_count; ++reaction) {
            const std::size_t row = balance_count + reaction;
            for (std::size_t supported = 0; supported < support_count; ++supported) {
                const std::size_t species = layout.supported_species[supported];
                const double coefficient = system.reaction_matrix(reaction, species);
                residual(static_cast<Eigen::Index>(row)) +=
                    coefficient * phase.chemical_potentials[species];
                for (std::size_t column = 0; column < dimension; ++column) {
                    for (std::size_t physical = 0;
                         physical < physical_dimension;
                         ++physical) {
                        const std::size_t full_physical = physical == support_count
                            ? system.species_count
                            : layout.supported_species[physical];
                        jacobian(
                            static_cast<Eigen::Index>(row),
                            static_cast<Eigen::Index>(column)
                        ) += coefficient
                            * phase.physical_hessian[
                                supported * physical_dimension + physical
                            ]
                            * phase.state_jacobian[
                                full_physical * dimension + column
                            ];
                    }
                }
            }
        }
        const std::size_t pressure_row = dimension - 1;
        residual(static_cast<Eigen::Index>(pressure_row)) =
            (phase.pressure - pressure_pa) / pressure_pa;
        for (std::size_t column = 0; column < dimension; ++column) {
            for (std::size_t physical = 0;
                 physical < physical_dimension;
                 ++physical) {
                const std::size_t full_physical = physical == support_count
                    ? system.species_count
                    : layout.supported_species[physical];
                jacobian(
                    static_cast<Eigen::Index>(pressure_row),
                    static_cast<Eigen::Index>(column)
                ) -= kGasConstantJPerMolK * temperature_k / pressure_pa
                    * phase.physical_hessian[
                        support_count * physical_dimension + physical
                    ]
                    * phase.state_jacobian[
                        full_physical * dimension + column
                    ];
            }
        }
        return std::pair{std::move(residual), std::move(jacobian)};
    };

    for (int iteration = 0; iteration < 10; ++iteration) {
        const auto [residual, jacobian] = residual_and_jacobian(variables);
        const double norm = residual.lpNorm<Eigen::Infinity>();
        if (!std::isfinite(norm)
            || norm <= 0.01 * kReactiveFinitePhaseSolverTolerance) {
            return;
        }
        Eigen::FullPivLU<Eigen::MatrixXd> factor(jacobian);
        if (static_cast<std::size_t>(factor.rank()) != dimension) {
            return;
        }
        const Eigen::VectorXd step = factor.solve(-residual);
        if (!step.allFinite()) {
            return;
        }
        bool accepted = false;
        for (int backtrack = 0; backtrack < 12; ++backtrack) {
            const double scale = std::ldexp(1.0, -backtrack);
            std::vector<double> candidate = variables;
            bool interior = true;
            for (std::size_t column = 0; column < dimension; ++column) {
                candidate[column] += scale * step(static_cast<Eigen::Index>(column));
                interior = interior
                    && candidate[column] > lower[column]
                    && candidate[column] < upper[column];
            }
            if (!interior) {
                continue;
            }
            try {
                const auto [candidate_residual, unused] =
                    residual_and_jacobian(candidate);
                if (candidate_residual.lpNorm<Eigen::Infinity>() < norm) {
                    variables = std::move(candidate);
                    accepted = true;
                    break;
                }
            } catch (const std::exception&) {
                continue;
            }
        }
        if (!accepted) {
            return;
        }
    }
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
    const std::vector<PhaseLayout>& layouts
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
        layouts
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
            layouts
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
            layouts
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
        std::vector<PhaseLayout> layouts
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
          layouts_(std::move(layouts)),
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
                layouts_
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
    std::vector<PhaseLayout> layouts_;
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
        return 0.0;
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
                    std::abs(
                        phases[left].amounts[species] / left_total
                        - phases[right].amounts[species] / right_total
                    )
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

std::pair<double, double> reduced_hessian_minimum(
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
        return {std::numeric_limits<double>::infinity(), 1.0};
    }
    const Eigen::MatrixXd nullspace = svd.matrixV().rightCols(
        static_cast<Eigen::Index>(variable_count - rank)
    );
    const Eigen::MatrixXd reduced = nullspace.transpose()
        * (0.5 * (hessian + hessian.transpose())) * nullspace;
    Eigen::VectorXd diagonal_scale(reduced.rows());
    for (Eigen::Index index = 0; index < reduced.rows(); ++index) {
        diagonal_scale(index) = 1.0 / std::sqrt(std::max(
            std::abs(reduced(index, index)),
            std::numeric_limits<double>::min()
        ));
    }
    const Eigen::MatrixXd equilibrated = diagonal_scale.asDiagonal()
        * reduced * diagonal_scale.asDiagonal();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(equilibrated);
    if (eigen.info() != Eigen::Success) {
        throw std::runtime_error("reactive finite-phase reduced Hessian eigensolve failed");
    }
    double scale = 1.0;
    for (Eigen::Index row = 0; row < equilibrated.rows(); ++row) {
        scale = std::max(scale, equilibrated.row(row).cwiseAbs().sum());
    }
    return {eigen.eigenvalues().minCoeff(), scale};
}

void populate_sensitivities(
    ReactiveFinitePhaseResult& result,
    const MultiphaseEvaluation& evaluation,
    const ProviderActiveParameterSet& active_parameters,
    const NormalizedBalances& balances,
    const std::vector<double>& g_ref_pressure_derivatives_per_pa,
    const std::vector<double>& g_ref_parameter_derivatives,
    const std::vector<double>& multipliers,
    double temperature_k,
    std::size_t species_count,
    const std::vector<PhaseLayout>& layouts
) {
    const std::size_t active_parameter_count = active_parameters.identities.size();
    const std::size_t phase_count = evaluation.phases.size();
    const std::size_t variable_count = evaluation.gradient.size();
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
    result.kkt_dimension = kkt_dimension;
    result.kkt_root_jacobian.assign(kkt_dimension * kkt_dimension, 0.0);
    for (std::size_t row = 0; row < kkt_dimension; ++row) {
        for (std::size_t column = 0; column < kkt_dimension; ++column) {
            result.kkt_root_jacobian[row * kkt_dimension + column] = kkt(
                static_cast<Eigen::Index>(row),
                static_cast<Eigen::Index>(column)
            );
        }
    }
    result.kkt_root_values = lagrangian_gradient(evaluation, multipliers);
    result.kkt_root_values.insert(
        result.kkt_root_values.end(),
        evaluation.constraints.begin(),
        evaluation.constraints.end()
    );
    Eigen::VectorXd row_scale = Eigen::VectorXd::Ones(
        static_cast<Eigen::Index>(kkt_dimension)
    );
    Eigen::VectorXd column_scale = Eigen::VectorXd::Ones(
        static_cast<Eigen::Index>(kkt_dimension)
    );
    for (int pass = 0; pass < 4; ++pass) {
        for (std::size_t row = 0; row < kkt_dimension; ++row) {
            double scale = 0.0;
            for (std::size_t column = 0; column < kkt_dimension; ++column) {
                scale = std::max(
                    scale,
                    std::abs(
                        kkt(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column))
                        * row_scale(static_cast<Eigen::Index>(row))
                        * column_scale(static_cast<Eigen::Index>(column))
                    )
                );
            }
            if (scale > 0.0 && std::isfinite(scale)) {
                row_scale(static_cast<Eigen::Index>(row)) /= scale;
            }
        }
        for (std::size_t column = 0; column < kkt_dimension; ++column) {
            double scale = 0.0;
            for (std::size_t row = 0; row < kkt_dimension; ++row) {
                scale = std::max(
                    scale,
                    std::abs(
                        kkt(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column))
                        * row_scale(static_cast<Eigen::Index>(row))
                        * column_scale(static_cast<Eigen::Index>(column))
                    )
                );
            }
            if (scale > 0.0 && std::isfinite(scale)) {
                column_scale(static_cast<Eigen::Index>(column)) /= scale;
            }
        }
    }
    const Eigen::MatrixXd scaled_kkt = row_scale.asDiagonal()
        * kkt * column_scale.asDiagonal();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(scaled_kkt);
    const double largest = svd.singularValues().size() == 0
        ? 0.0
        : svd.singularValues()(0);
    const double threshold = 256.0 * std::numeric_limits<double>::epsilon()
        * static_cast<double>(std::max<std::size_t>(1, kkt_dimension))
        * std::max(1.0, largest);
    result.kkt_rank = 0;
    double smallest = std::numeric_limits<double>::infinity();
    for (Eigen::Index index = 0; index < svd.singularValues().size(); ++index) {
        const double value = svd.singularValues()(index);
        if (value > threshold) {
            ++result.kkt_rank;
            smallest = std::min(smallest, value);
        }
    }
    result.condition_number_inf = result.kkt_rank == kkt_dimension && smallest > 0.0
        ? largest / smallest
        : std::numeric_limits<double>::infinity();
    Eigen::FullPivLU<Eigen::MatrixXd> factor(scaled_kkt);
    if (result.kkt_rank != kkt_dimension) {
        result.sensitivity_status = "unavailable";
        result.sensitivity_failure_reason = "singular_exact_kkt_system";
        return;
    }
    if (!std::isfinite(result.condition_number_inf)
        || result.condition_number_inf > kChemicalConditionNumberMax) {
        result.sensitivity_status = "unavailable";
        result.sensitivity_failure_reason = "ill_conditioned_exact_kkt_system";
        return;
    }
    if (!g_ref_pressure_derivatives_per_pa.empty()
        && g_ref_pressure_derivatives_per_pa.size() != species_count) {
        result.sensitivity_status = "unavailable";
        result.sensitivity_failure_reason = "incomplete_reference_pressure_derivatives";
        return;
    }
    const std::size_t input_count = constraint_count + 1 + active_parameter_count;
    result.input_order.reserve(input_count);
    for (std::size_t row = 0; row < constraint_count; ++row) {
        result.input_order.push_back("balance_total[" + std::to_string(row) + "]");
    }
    result.input_order.emplace_back("pressure_pa");
    result.input_order.insert(
        result.input_order.end(),
        active_parameters.identities.begin(),
        active_parameters.identities.end()
    );
    result.kkt_input_jacobian.assign(kkt_dimension * input_count, 0.0);
    result.input_amount_derivatives.assign(
        phase_count * species_count * input_count, 0.0
    );
    result.input_volume_derivatives.assign(phase_count * input_count, 0.0);
    result.input_pressure_derivatives.assign(phase_count * input_count, 0.0);
    result.input_chemical_potential_derivatives.assign(
        phase_count * species_count * input_count, 0.0
    );
    result.input_log_fugacity_derivatives.assign(
        phase_count * species_count * input_count, 0.0
    );
    for (std::size_t input = 0; input < input_count; ++input) {
        Eigen::VectorXd partial = Eigen::VectorXd::Zero(
            static_cast<Eigen::Index>(kkt_dimension)
        );
        if (input < constraint_count) {
            partial(static_cast<Eigen::Index>(variable_count + input)) =
                -1.0 / balances.scales[input];
        } else {
            const bool pressure_input = input == constraint_count;
            const std::size_t active = pressure_input
                ? 0
                : input - constraint_count - 1;
            for (std::size_t phase = 0; phase < phase_count; ++phase) {
                const PhaseLayout& layout = layouts[phase];
                for (std::size_t coordinate = 0;
                     coordinate < layout.dimension();
                     ++coordinate) {
                    partial(static_cast<Eigen::Index>(layout.offset + coordinate)) =
                        pressure_input
                        ? evaluation.phases[phase].objective_pressure_cross[coordinate]
                        : evaluation.phases[phase].objective_parameter_cross[
                            coordinate * active_parameter_count + active
                        ];
                }
            }
        }
        for (std::size_t row = 0; row < kkt_dimension; ++row) {
            result.kkt_input_jacobian[row * input_count + input] =
                partial(static_cast<Eigen::Index>(row));
        }
        const Eigen::VectorXd response = column_scale.asDiagonal()
            * factor.solve(row_scale.asDiagonal() * -partial);
        if (!response.allFinite()) {
            result.sensitivity_status = "unavailable";
            result.sensitivity_failure_reason = "nonfinite_exact_kkt_response";
            return;
        }
        for (std::size_t phase = 0; phase < phase_count; ++phase) {
            const PhaseBlock& block = evaluation.phases[phase];
            const PhaseLayout& layout = layouts[phase];
            const std::size_t offset = layout.offset;
            const std::size_t phase_dimension = layout.dimension();
            const std::size_t support_count = layout.supported_species.size();
            const std::size_t physical_dimension = support_count + 1;
            std::vector<double> full_physical_response(species_count + 1, 0.0);
            for (std::size_t physical = 0;
                 physical < species_count + 1;
                 ++physical) {
                for (std::size_t coordinate = 0; coordinate < phase_dimension; ++coordinate) {
                    full_physical_response[physical] +=
                        block.state_jacobian[physical * phase_dimension + coordinate]
                        * response(static_cast<Eigen::Index>(offset + coordinate));
                }
            }
            const bool pressure_input = input == constraint_count;
            const bool active_input = input > constraint_count;
            const std::size_t active = active_input
                ? input - constraint_count - 1
                : 0;
            if (active_input) {
                full_physical_response.back() += block.volume_parameter_direct[active];
            }
            for (std::size_t species = 0; species < species_count; ++species) {
                result.input_amount_derivatives[
                    (phase * species_count + species) * input_count + input
                ] = full_physical_response[species];
            }
            for (std::size_t supported = 0; supported < support_count; ++supported) {
                const std::size_t species = layout.supported_species[supported];
                double reference_direct = 0.0;
                double chemical_response = 0.0;
                if (pressure_input && !g_ref_pressure_derivatives_per_pa.empty()) {
                    reference_direct = g_ref_pressure_derivatives_per_pa[species];
                    chemical_response = reference_direct;
                } else if (active_input) {
                    chemical_response = block.chemical_potential_parameter_direct[
                        species * active_parameter_count + active
                    ];
                    reference_direct = g_ref_parameter_derivatives.empty()
                        ? 0.0
                        : g_ref_parameter_derivatives[
                            species * active_parameter_count + active
                        ];
                }
                for (std::size_t physical = 0;
                     physical < physical_dimension;
                     ++physical) {
                    chemical_response +=
                        block.physical_hessian[
                            supported * physical_dimension + physical
                        ]
                        * (physical == support_count
                            ? full_physical_response.back()
                            : full_physical_response[
                                layout.supported_species[physical]
                            ]);
                }
                result.input_chemical_potential_derivatives[
                    (phase * species_count + species) * input_count + input
                ] = chemical_response;
                result.input_log_fugacity_derivatives[
                    (phase * species_count + species) * input_count + input
                ] = chemical_response - reference_direct;
            }
            result.input_volume_derivatives[phase * input_count + input] =
                full_physical_response.back();
            double pressure_response = active_input
                ? block.pressure_parameter_direct[active]
                : 0.0;
            for (std::size_t physical = 0;
                 physical < physical_dimension;
                 ++physical) {
                pressure_response -= kGasConstantJPerMolK * temperature_k
                    * block.physical_hessian[
                        support_count * physical_dimension + physical
                    ] * (physical == support_count
                        ? full_physical_response.back()
                        : full_physical_response[layout.supported_species[physical]]);
            }
            result.input_pressure_derivatives[phase * input_count + input] =
                pressure_response;
        }
    }
    result.sensitivity_status = "available";
    result.sensitivity_failure_reason.clear();
    result.parameter_order = active_parameters.identities;
    result.amount_derivatives.assign(
        phase_count * species_count * active_parameter_count, 0.0
    );
    result.volume_derivatives.assign(phase_count * active_parameter_count, 0.0);
    result.pressure_derivatives.assign(phase_count * active_parameter_count, 0.0);
    result.chemical_potential_derivatives.assign(
        phase_count * species_count * active_parameter_count, 0.0
    );
    result.log_fugacity_derivatives.assign(
        phase_count * species_count * active_parameter_count, 0.0
    );
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        for (std::size_t active = 0; active < active_parameter_count; ++active) {
            const std::size_t input = constraint_count + 1 + active;
            result.volume_derivatives[phase * active_parameter_count + active] =
                result.input_volume_derivatives[phase * input_count + input];
            result.pressure_derivatives[phase * active_parameter_count + active] =
                result.input_pressure_derivatives[phase * input_count + input];
            for (std::size_t species = 0; species < species_count; ++species) {
                const std::size_t source =
                    (phase * species_count + species) * input_count + input;
                const std::size_t target =
                    (phase * species_count + species) * active_parameter_count + active;
                result.amount_derivatives[target] = result.input_amount_derivatives[source];
                result.chemical_potential_derivatives[target] =
                    result.input_chemical_potential_derivatives[source];
                result.log_fugacity_derivatives[target] =
                    result.input_log_fugacity_derivatives[source];
            }
        }
    }
}

std::pair<double, std::vector<double>> contract_output(
    const ReactiveFinitePhaseOutput& output,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& amounts,
    double volume,
    const std::vector<double>& chemical_potentials,
    const std::vector<double>& log_fugacities,
    const std::vector<double>& amount_derivatives,
    const std::vector<double>& volume_derivatives,
    const std::vector<double>& pressure_derivatives,
    const std::vector<double>& chemical_potential_derivatives,
    const std::vector<double>& log_fugacity_derivatives,
    std::size_t derivative_count
) {
    const std::size_t species_count = amounts.size();
    if (output.coefficients.size() != species_count
        || (!output.solvent_mass_coefficients_kg_per_mol.empty()
            && output.solvent_mass_coefficients_kg_per_mol.size() != species_count)
        || amount_derivatives.size() != species_count * derivative_count
        || volume_derivatives.size() != derivative_count
        || pressure_derivatives.size() != derivative_count
        || chemical_potential_derivatives.size() != species_count * derivative_count
        || log_fugacity_derivatives.size() != species_count * derivative_count) {
        throw std::invalid_argument("reactive finite-phase output dimensions changed");
    }
    std::vector<double> derivatives(derivative_count, 0.0);
    if (output.selector == "system.pressure" || output.selector == "phase.pressure") {
        derivatives = pressure_derivatives;
        return {pressure_pa, std::move(derivatives)};
    }
    const double total = std::accumulate(amounts.begin(), amounts.end(), 0.0);
    const double numerator = std::inner_product(
        output.coefficients.begin(), output.coefficients.end(), amounts.begin(), 0.0
    );
    if (output.selector == "phase.amount") {
        for (std::size_t derivative = 0; derivative < derivative_count; ++derivative) {
            for (std::size_t species = 0; species < species_count; ++species) {
                derivatives[derivative] += output.coefficients[species]
                    * amount_derivatives[species * derivative_count + derivative];
            }
        }
        return {numerator, std::move(derivatives)};
    }
    if (output.selector == "phase.mole_fraction") {
        for (std::size_t derivative = 0; derivative < derivative_count; ++derivative) {
            double dnumerator = 0.0;
            double dtotal = 0.0;
            for (std::size_t species = 0; species < species_count; ++species) {
                const double value = amount_derivatives[
                    species * derivative_count + derivative
                ];
                dnumerator += output.coefficients[species] * value;
                dtotal += value;
            }
            derivatives[derivative] =
                (dnumerator * total - numerator * dtotal) / (total * total);
        }
        return {numerator / total, std::move(derivatives)};
    }
    if (output.selector == "phase.molality") {
        const double denominator = std::inner_product(
            output.solvent_mass_coefficients_kg_per_mol.begin(),
            output.solvent_mass_coefficients_kg_per_mol.end(),
            amounts.begin(),
            0.0
        );
        if (denominator <= 0.0) {
            throw std::domain_error("reactive finite-phase molality denominator is not positive");
        }
        for (std::size_t derivative = 0; derivative < derivative_count; ++derivative) {
            double dnumerator = 0.0;
            double ddenominator = 0.0;
            for (std::size_t species = 0; species < species_count; ++species) {
                const double value = amount_derivatives[
                    species * derivative_count + derivative
                ];
                dnumerator += output.coefficients[species] * value;
                ddenominator += output.solvent_mass_coefficients_kg_per_mol[species]
                    * value;
            }
            derivatives[derivative] =
                (dnumerator * denominator - numerator * ddenominator)
                / (denominator * denominator);
        }
        return {numerator / denominator, std::move(derivatives)};
    }
    if (output.selector == "phase.molar_density") {
        for (std::size_t derivative = 0; derivative < derivative_count; ++derivative) {
            double dtotal = 0.0;
            for (std::size_t species = 0; species < species_count; ++species) {
                dtotal += amount_derivatives[species * derivative_count + derivative];
            }
            derivatives[derivative] = dtotal / volume
                - total * volume_derivatives[derivative] / (volume * volume);
        }
        return {total / volume, std::move(derivatives)};
    }
    if (output.selector == "phase.compressibility_factor") {
        const double value = pressure_pa * volume
            / (total * kGasConstantJPerMolK * temperature_k);
        for (std::size_t derivative = 0; derivative < derivative_count; ++derivative) {
            double dtotal = 0.0;
            for (std::size_t species = 0; species < species_count; ++species) {
                dtotal += amount_derivatives[species * derivative_count + derivative];
            }
            derivatives[derivative] = value * (
                pressure_derivatives[derivative] / pressure_pa
                + volume_derivatives[derivative] / volume
                - dtotal / total
            );
        }
        return {value, std::move(derivatives)};
    }
    const bool fugacity = output.selector == "phase.fugacity";
    const bool log_fugacity = output.selector == "phase.log_fugacity" || fugacity;
    if (!log_fugacity && output.selector != "phase.chemical_potential_over_rt") {
        throw std::invalid_argument("unsupported reactive finite-phase output selector");
    }
    const std::vector<double>& values = log_fugacity
        ? log_fugacities
        : chemical_potentials;
    const std::vector<double>& source_derivatives = log_fugacity
        ? log_fugacity_derivatives
        : chemical_potential_derivatives;
    double value = std::inner_product(
        output.coefficients.begin(),
        output.coefficients.end(),
        values.begin(),
        0.0,
        std::plus<>(),
        [](double coefficient, double component) {
            return coefficient == 0.0 ? 0.0 : coefficient * component;
        }
    );
    for (std::size_t derivative = 0; derivative < derivative_count; ++derivative) {
        for (std::size_t species = 0; species < species_count; ++species) {
            if (output.coefficients[species] != 0.0) {
                derivatives[derivative] += output.coefficients[species]
                    * source_derivatives[species * derivative_count + derivative];
            }
        }
    }
    if (fugacity) {
        value = std::exp(value);
        for (double& derivative : derivatives) {
            derivative *= value;
        }
    }
    return {value, std::move(derivatives)};
}

void populate_outputs(
    ReactiveFinitePhaseResult& result,
    const MultiphaseEvaluation& evaluation,
    const ProviderActiveParameterSet& active_parameters,
    const std::vector<double>& g_ref_parameter_derivatives,
    double temperature_k,
    double pressure_pa,
    double amount_scale,
    const std::vector<PhaseLayout>& layouts,
    const std::vector<ReactiveFinitePhaseOutput>& outputs
) {
    const std::size_t phase_count = evaluation.phases.size();
    const std::size_t species_count = evaluation.phases.front().amounts.size();
    const std::size_t input_count = result.input_order.size();
    const std::size_t active_count = active_parameters.identities.size();
    const std::size_t primal_count = evaluation.gradient.size();
    const std::size_t local_count = primal_count + evaluation.constraints.size();
    result.output_values.reserve(outputs.size());
    result.output_input_derivatives.reserve(outputs.size() * input_count);
    result.output_local_derivatives.reserve(outputs.size() * local_count);
    result.output_parameter_partials.reserve(outputs.size() * active_count);
    for (const ReactiveFinitePhaseOutput& output : outputs) {
        if (output.phase_index >= phase_count) {
            throw std::invalid_argument("reactive finite-phase output names an unknown phase");
        }
        const std::size_t phase = output.phase_index;
        const PhaseBlock& block = evaluation.phases[phase];
        const PhaseLayout& layout = layouts[phase];
        const std::vector<double>& amounts = result.phase_amounts_mol[phase];
        const double volume = result.phase_volumes_m3[phase];
        std::vector<double> input_amounts(species_count * input_count, 0.0);
        std::vector<double> input_mu(species_count * input_count, 0.0);
        std::vector<double> input_log_f(species_count * input_count, 0.0);
        for (std::size_t species = 0; species < species_count; ++species) {
            for (std::size_t input = 0; input < input_count; ++input) {
                const std::size_t source =
                    (phase * species_count + species) * input_count + input;
                input_amounts[species * input_count + input] =
                    result.input_amount_derivatives[source];
                input_mu[species * input_count + input] =
                    result.input_chemical_potential_derivatives[source];
                input_log_f[species * input_count + input] =
                    result.input_log_fugacity_derivatives[source];
            }
        }
        std::vector<double> input_volumes(input_count, 0.0);
        std::vector<double> input_pressures(input_count, 0.0);
        for (std::size_t input = 0; input < input_count; ++input) {
            input_volumes[input] =
                result.input_volume_derivatives[phase * input_count + input];
            input_pressures[input] = result.input_order[input] == "pressure_pa"
                ? 1.0
                : 0.0;
        }
        const auto input = contract_output(
            output,
            temperature_k,
            pressure_pa,
            amounts,
            volume,
            result.phase_chemical_potentials_over_rt[phase],
            result.phase_log_fugacities[phase],
            input_amounts,
            input_volumes,
            input_pressures,
            input_mu,
            input_log_f,
            input_count
        );
        result.output_values.push_back(input.first);
        result.output_input_derivatives.insert(
            result.output_input_derivatives.end(), input.second.begin(), input.second.end()
        );

        std::vector<double> local_amounts(species_count * local_count, 0.0);
        std::vector<double> local_volumes(local_count, 0.0);
        std::vector<double> local_pressures(local_count, 0.0);
        std::vector<double> local_mu(species_count * local_count, 0.0);
        std::vector<double> local_log_f(species_count * local_count, 0.0);
        const std::size_t support_count = layout.supported_species.size();
        const std::size_t physical_dimension = support_count + 1;
        for (std::size_t coordinate = 0; coordinate < layout.dimension(); ++coordinate) {
            const std::size_t local = layout.offset + coordinate;
            std::vector<double> physical_response(species_count + 1, 0.0);
            for (std::size_t physical = 0; physical < species_count + 1; ++physical) {
                physical_response[physical] = amount_scale
                    * block.state_jacobian[physical * layout.dimension() + coordinate];
            }
            for (std::size_t species = 0; species < species_count; ++species) {
                local_amounts[species * local_count + local] = physical_response[species];
            }
            local_volumes[local] = physical_response.back();
            for (std::size_t supported = 0; supported < support_count; ++supported) {
                const std::size_t species = layout.supported_species[supported];
                double derivative = 0.0;
                for (std::size_t physical = 0; physical < physical_dimension; ++physical) {
                    derivative += block.physical_hessian[
                        supported * physical_dimension + physical
                    ] * (physical == support_count
                        ? physical_response.back() / amount_scale
                        : physical_response[layout.supported_species[physical]] / amount_scale);
                }
                local_mu[species * local_count + local] = derivative;
                local_log_f[species * local_count + local] = derivative;
            }
            double derivative = 0.0;
            for (std::size_t physical = 0; physical < physical_dimension; ++physical) {
                derivative -= kGasConstantJPerMolK * temperature_k
                    * block.physical_hessian[
                        support_count * physical_dimension + physical
                    ] * (physical == support_count
                        ? physical_response.back() / amount_scale
                        : physical_response[layout.supported_species[physical]] / amount_scale);
            }
            local_pressures[local] = derivative;
        }
        const auto local = contract_output(
            output,
            temperature_k,
            pressure_pa,
            amounts,
            volume,
            result.phase_chemical_potentials_over_rt[phase],
            result.phase_log_fugacities[phase],
            local_amounts,
            local_volumes,
            local_pressures,
            local_mu,
            local_log_f,
            local_count
        );
        result.output_local_derivatives.insert(
            result.output_local_derivatives.end(), local.second.begin(), local.second.end()
        );

        std::vector<double> parameter_amounts(species_count * active_count, 0.0);
        std::vector<double> parameter_volumes(active_count, 0.0);
        std::vector<double> parameter_pressures(active_count, 0.0);
        std::vector<double> parameter_mu(species_count * active_count, 0.0);
        std::vector<double> parameter_log_f(species_count * active_count, 0.0);
        for (std::size_t active = 0; active < active_count; ++active) {
            const double volume_direct = block.volume_parameter_direct[active];
            parameter_volumes[active] = amount_scale * volume_direct;
            double pressure_direct = block.pressure_parameter_direct[active];
            for (std::size_t supported = 0; supported < support_count; ++supported) {
                const std::size_t species = layout.supported_species[supported];
                const double mu_direct = block.chemical_potential_parameter_direct[
                    species * active_count + active
                ] + block.physical_hessian[
                    supported * physical_dimension + support_count
                ] * volume_direct;
                parameter_mu[species * active_count + active] = mu_direct;
                parameter_log_f[species * active_count + active] = mu_direct
                    - (g_ref_parameter_derivatives.empty()
                        ? 0.0
                        : g_ref_parameter_derivatives[species * active_count + active]);
            }
            pressure_direct -= kGasConstantJPerMolK * temperature_k
                * block.physical_hessian.back() * volume_direct;
            parameter_pressures[active] = pressure_direct;
        }
        const auto parameter = contract_output(
            output,
            temperature_k,
            pressure_pa,
            amounts,
            volume,
            result.phase_chemical_potentials_over_rt[phase],
            result.phase_log_fugacities[phase],
            parameter_amounts,
            parameter_volumes,
            parameter_pressures,
            parameter_mu,
            parameter_log_f,
            active_count
        );
        result.output_parameter_partials.insert(
            result.output_parameter_partials.end(),
            parameter.second.begin(),
            parameter.second.end()
        );
    }
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
    const PhaseDirectionCompilation* phase_directions,
    const std::vector<std::vector<std::size_t>>& phase_support_indices,
    const std::vector<ReactiveFinitePhaseOutput>& outputs
) {
    ReactiveFinitePhaseResult result;
    const std::size_t phase_count = starts.size();
    const std::size_t species_count = system.species_count;
    result.resolved_policy = {
        {"ipopt.tol", kReactiveFinitePhaseSolverTolerance},
        {"ipopt.constr_viol_tol", kReactiveFinitePhaseSolverTolerance},
        {"ipopt.acceptable_tol", kReactiveFinitePhaseAcceptableTolerance},
        {"ipopt.max_iter", static_cast<double>(max_iterations)},
        {"start.balance_inf_norm_max", kReactiveFinitePhaseStartBalanceTolerance},
        {"certificate.balance_inf_norm_max", kChemicalBalanceTolerance},
        {"certificate.pressure_relative_inf_norm_max", kChemicalPressureTolerance},
        {"certificate.affinity_inf_norm_max", kChemicalAffinityTolerance},
        {"certificate.kkt_stationarity_inf_norm_max", kChemicalKktTolerance},
        {"certificate.minimum_amount_mol_min_exclusive", trace_floor},
        {"certificate.bound_margin_min_exclusive", kReactiveFinitePhaseBoundMarginMinimum},
        {"certificate.strict_minimum_scale_factor", kReactiveFinitePhaseStrictMinimumFactor},
    };
    if (phase_count > 1) {
        result.resolved_policy.push_back({
            "certificate.phase_distinctness_min_exclusive",
            kReactiveFinitePhaseDistinctnessMinimum,
        });
    }
    if (phase_count < 1 || packing_fraction_intervals.size() != phase_count
        || phase_support_indices.size() != phase_count) {
        throw std::invalid_argument(
            "reactive finite-phase solve requires at least one declared phase"
        );
    }
    const std::vector<PhaseLayout> layouts = make_phase_layouts(
        system.charges, phase_support_indices
    );
    const std::size_t variable_count =
        layouts.back().offset + layouts.back().dimension();
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
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        std::vector<bool> supported(species_count, false);
        for (std::size_t species : layouts[phase].supported_species) {
            supported[species] = true;
        }
        for (std::size_t species = 0; species < species_count; ++species) {
            if (supported[species]) {
                continue;
            }
            if (starts[phase].amounts_mol.size() == species_count
                && starts[phase].amounts_mol[species] != 0.0) {
                throw std::invalid_argument(
                    "reactive finite-phase start populates an unsupported species"
                );
            }
            if (phase_directions != nullptr) {
                for (std::size_t direction = 0;
                     direction < phase_directions->direction_matrix.rows;
                     ++direction) {
                    if (phase_directions->direction_matrix(
                            direction, phase * species_count + species
                        ) != 0.0) {
                        throw std::invalid_argument(
                            "phase direction populates an unsupported species"
                        );
                    }
                }
            } else {
                for (std::size_t reaction = 0;
                     reaction < system.reaction_matrix.rows;
                     ++reaction) {
                    if (system.reaction_matrix(reaction, species) != 0.0) {
                        throw std::invalid_argument(
                            "sparse reactive phases require declared phase directions"
                        );
                    }
                }
            }
        }
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
    if (system.ideal_gas && !parameters.identities.empty()) {
        throw std::invalid_argument(
            "reactive ideal-gas phases do not expose Provider active parameters"
        );
    }
    std::vector<double> initial(variable_count, 0.0);
    std::vector<double> lower(variable_count, 0.0);
    std::vector<double> upper(variable_count, 0.0);
    // The compiler inserts the positive molar-mass row for species reactions.
    // Configured phase directions are admitted only when they conserve the same
    // mass functional and their reference phases sum to the feed.  Therefore
    // m_i n_i <= sum_j m_j n_j = feed mass for every nonnegative admitted state;
    // twice that ceiling is a rigorous guard with the physical boundary inactive.
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
        const PhaseLayout& layout = layouts[phase];
        const AmountChart& chart = layout.chart;
        const std::size_t amount_dimension = layout.amount_dimension();
        std::vector<bool> supported(species_count, false);
        for (std::size_t species : layout.supported_species) {
            supported[species] = true;
        }
        std::vector<double> start_amounts = start.amounts_mol;
        bool valid_amounts = start_amounts.size() == species_count;
        for (std::size_t species = 0;
             valid_amounts && species < species_count;
             ++species) {
            const double amount = start_amounts[species];
            valid_amounts = std::isfinite(amount)
                && (supported[species] ? amount > trace_floor : amount == 0.0);
        }
        if (!valid_amounts && phase_count == 1) {
            const std::optional<std::vector<double>> seed = strict_positive_seed(
                system,
                layout.supported_species,
                trace_floor,
                provider.metadata().ion_mole_fraction_max
            );
            if (seed.has_value()) {
                start_amounts = *seed;
                valid_amounts = true;
            }
        }
        if (!valid_amounts
            || !std::isfinite(start.molar_volume_m3_per_mol)
            || start.molar_volume_m3_per_mol <= 0.0
            || (!system.ideal_gas
                && (!std::isfinite(packing[0]) || !std::isfinite(packing[1])
                    || packing[0] <= 0.0 || packing[1] <= packing[0]
                    || packing[1] >= 1.0))) {
            throw std::invalid_argument("reactive finite-phase start or domain is invalid");
        }
        const std::size_t offset = layout.offset;
        const double total = std::accumulate(
            start_amounts.begin(), start_amounts.end(), 0.0
        );
        const double volume = total * start.molar_volume_m3_per_mol;
        const double start_packing = system.ideal_gas
            ? 0.0
            : provider.evaluate_packing_fraction(
                temperature_k, start_amounts, volume
            ).value;
        if (!system.ideal_gas
            && (start_packing <= packing[0] || start_packing >= packing[1])) {
            throw std::invalid_argument(
                "reactive finite-phase start lies outside the declared packing interval"
            );
        }
        std::vector<double> supported_amounts;
        supported_amounts.reserve(layout.supported_species.size());
        for (std::size_t species : layout.supported_species) {
            supported_amounts.push_back(start_amounts[species]);
        }
        const std::vector<double> amount_coordinates =
            invert_amount_chart(chart, supported_amounts);
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
            for (std::size_t supported_index = 0;
                 supported_index < layout.supported_species.size();
                 ++supported_index) {
                const std::size_t species =
                    layout.supported_species[supported_index];
                lower[offset + supported_index] = std::log(0.1 * trace_floor);
                upper[offset + supported_index] = std::log(std::nextafter(
                    2.0 * total_feed_mass_kg
                        / system.molar_masses_kg_per_mol[species],
                    std::numeric_limits<double>::infinity()
                ));
            }
        } else {
            double charge_equivalent_upper = 0.0;
            for (std::size_t supported_index : chart.cation_indices) {
                const std::size_t species =
                    layout.supported_species[supported_index];
                charge_equivalent_upper +=
                    static_cast<double>(chart.charges[supported_index])
                    * 2.0 * total_feed_mass_kg
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
                const std::size_t species = layout.supported_species[
                    chart.neutral_indices[neutral]
                ];
                lower[offset + neutral_offset + neutral] =
                    std::log(0.1 * trace_floor);
                    upper[offset + neutral_offset + neutral] =
                        std::log(std::nextafter(
                        2.0 * total_feed_mass_kg
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
        if (system.ideal_gas) {
            const double minimum_volume = 0.01 * trace_floor
                * kGasConstantJPerMolK * temperature_k
                / (pressure_pa * total_feed_amount_mol);
            const double maximum_amount = total_feed_mass_kg
                / *std::min_element(
                    system.molar_masses_kg_per_mol.begin(),
                    system.molar_masses_kg_per_mol.end()
                );
            const double maximum_volume = 100.0 * maximum_amount
                * kGasConstantJPerMolK * temperature_k
                / (pressure_pa * total_feed_amount_mol);
            initial[offset + amount_dimension] = std::log(
                volume / total_feed_amount_mol
            );
            lower[offset + amount_dimension] = std::log(minimum_volume);
            upper[offset + amount_dimension] = std::log(maximum_volume);
        } else {
            initial[offset + amount_dimension] = std::log(start_packing);
            lower[offset + amount_dimension] = std::log(packing[0]);
            upper[offset + amount_dimension] = std::nextafter(
                std::log(packing[1]), -std::numeric_limits<double>::infinity()
            );
        }
    }

    NormalizedBalances balances = normalized_balances(
        system, phase_directions, phase_support_indices
    );
    for (double& coefficient : balances.matrix.values) {
        coefficient *= total_feed_amount_mol;
    }
    const std::size_t constraint_count = balances.matrix.rows;
    const std::size_t amount_variable_count = std::accumulate(
        layouts.begin(), layouts.end(), std::size_t{0},
        [](std::size_t total, const PhaseLayout& layout) {
            return total + layout.amount_dimension();
        }
    );
    if (amount_variable_count < constraint_count) {
        throw std::invalid_argument(
            "configured reactive phases have more balances than amount coordinates"
        );
    }
    result.admissible_nullity = amount_variable_count - constraint_count;
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
        layouts
    );
    if (inf_norm(initial_evaluation.constraints)
        > kReactiveFinitePhaseStartBalanceTolerance) {
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
        layouts
    );
    Ipopt::SmartPtr<Ipopt::TNLP> problem = raw_problem;
    Ipopt::SmartPtr<Ipopt::IpoptApplication> application = IpoptApplicationFactory();
    application->Options()->SetIntegerValue("print_level", 0);
    application->Options()->SetStringValue("sb", "yes");
    application->Options()->SetStringValue("hessian_approximation", "exact");
    application->Options()->SetStringValue("mu_strategy", "adaptive");
    application->Options()->SetStringValue("nlp_scaling_method", "none");
    application->Options()->SetNumericValue("tol", kReactiveFinitePhaseSolverTolerance);
    application->Options()->SetNumericValue(
        "constr_viol_tol", kReactiveFinitePhaseSolverTolerance
    );
    application->Options()->SetNumericValue(
        "acceptable_tol", kReactiveFinitePhaseAcceptableTolerance
    );
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
    if (phase_directions == nullptr) {
        polish_single_phase_reaction_root(
            provider,
            parameters,
            system,
            balances,
            temperature_k,
            pressure_pa,
            packing_fraction_intervals,
            layouts,
            lower,
            upper,
            solution
        );
    }
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
            layouts
        );
        double maximum_relative_residual = 0.0;
        bool changed = false;
        for (std::size_t phase = 0; phase < phase_count; ++phase) {
            const PhaseBlock& block = pressure_evaluation.phases[phase];
            const PhaseLayout& layout = layouts[phase];
            const std::size_t phase_dimension = layout.dimension();
            const std::size_t amount_dimension = layout.amount_dimension();
            const std::size_t support_count = layout.supported_species.size();
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
                    support_count * (support_count + 1) + support_count
                ]
                * volume_derivative;
            if (!std::isfinite(pressure_derivative)
                || std::abs(pressure_derivative)
                    <= std::numeric_limits<double>::min()) {
                continue;
            }
            const std::size_t coordinate = layout.offset + amount_dimension;
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
        layouts
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
        layouts
    );
    result.objective_over_rt = total_feed_amount_mol * evaluated.objective;
    result.balance_inf_norm = inf_norm(evaluated.constraints);
    result.phase_distinctness = phase_distinctness(evaluated.phases);
    result.phase_amounts_mol.reserve(phase_count);
    result.phase_volumes_m3.reserve(phase_count);
    result.phase_pressures_pa.reserve(phase_count);
    result.phase_packing_fractions.reserve(phase_count);
    result.phase_chemical_potentials_over_rt.reserve(phase_count);
    result.phase_log_fugacities.reserve(phase_count);
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
        const PhaseLayout& layout = layouts[phase];
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
        std::vector<double> log_fugacities = block.chemical_potentials;
        for (std::size_t species = 0; species < species_count; ++species) {
            log_fugacities[species] += std::log(
                kGasConstantJPerMolK * temperature_k
            ) - system.g_ref[species];
        }
        result.phase_log_fugacities.push_back(std::move(log_fugacities));
        result.pressure_relative_inf_norm = std::max(
            result.pressure_relative_inf_norm,
            std::abs(block.pressure - pressure_pa) / pressure_pa
        );
        for (std::size_t species : layout.supported_species) {
            minimum_amount = std::min(
                minimum_amount, total_feed_amount_mol * block.amounts[species]
            );
        }
        const std::size_t offset = layout.offset;
        const std::size_t phase_dimension = layout.dimension();
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
                    const double coefficient =
                        system.reaction_matrix(reaction, species);
                    if (coefficient != 0.0) {
                        affinity += coefficient * block.chemical_potentials[species];
                    }
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
                    const double coefficient = phase_directions->direction_matrix(
                        direction, phase * species_count + species
                    );
                    if (coefficient != 0.0) {
                        affinity += coefficient
                            * evaluated.phases[phase].chemical_potentials[species];
                    }
                }
            }
            result.reaction_affinity_inf_norm = std::max(
                result.reaction_affinity_inf_norm, std::abs(affinity)
            );
        }
        result.transfer_chemical_potential_inf_norm =
            result.reaction_affinity_inf_norm;
    }
    for (std::size_t phase = 1;
         phase_directions == nullptr && phase < phase_count;
         ++phase) {
        double electric_offset = 0.0;
        double charge_norm_squared = 0.0;
        std::vector<bool> shared(species_count, false);
        for (std::size_t species : layouts[0].supported_species) {
            shared[species] = std::binary_search(
                layouts[phase].supported_species.begin(),
                layouts[phase].supported_species.end(),
                species
            );
            if (shared[species]) {
                charge_norm_squared += static_cast<double>(system.charges[species])
                    * static_cast<double>(system.charges[species]);
            }
        }
        if (charge_norm_squared > 0.0) {
            for (std::size_t species = 0; species < species_count; ++species) {
                if (shared[species]) {
                    electric_offset += static_cast<double>(system.charges[species])
                    * (evaluated.phases[phase].chemical_potentials[species]
                        - evaluated.phases[0].chemical_potentials[species]);
                }
            }
            electric_offset /= charge_norm_squared;
        }
        for (std::size_t species = 0; species < species_count; ++species) {
            if (!shared[species]) {
                continue;
            }
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
    const auto [minimum_reduced_eigenvalue, reduced_hessian_scale] =
        reduced_hessian_minimum(
        evaluated, variable_count, evaluated.constraints.size()
    );
    result.reduced_hessian_min_eigenvalue = minimum_reduced_eigenvalue;
    const double strict_minimum_threshold =
        kReactiveFinitePhaseStrictMinimumFactor * reduced_hessian_scale;
    result.minimum_amount_mol = minimum_amount;
    result.minimum_bound_margin = minimum_bound_margin;
    result.reduced_hessian_scale = reduced_hessian_scale;
    result.reduced_hessian_threshold = strict_minimum_threshold;
    const bool physical_certificates_passed =
        result.balance_inf_norm <= kChemicalBalanceTolerance
        && result.pressure_relative_inf_norm <= kChemicalPressureTolerance
        && result.reaction_affinity_inf_norm <= kChemicalAffinityTolerance
        && result.transfer_chemical_potential_inf_norm <= kChemicalAffinityTolerance
        && result.kkt_stationarity_inf_norm <= kChemicalKktTolerance
        && (phase_count == 1
            || result.phase_distinctness > kReactiveFinitePhaseDistinctnessMinimum)
        && minimum_amount > trace_floor
        && minimum_bound_margin > kReactiveFinitePhaseBoundMarginMinimum
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
                layouts
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
        balances,
        system.g_ref_pressure_derivatives_per_pa,
        system.g_ref_parameter_derivatives,
        multipliers,
        temperature_k,
        species_count,
        layouts
    );
    for (double& derivative : result.amount_derivatives) {
        derivative *= total_feed_amount_mol;
    }
    for (double& derivative : result.volume_derivatives) {
        derivative *= total_feed_amount_mol;
    }
    for (double& derivative : result.input_amount_derivatives) {
        derivative *= total_feed_amount_mol;
    }
    for (double& derivative : result.input_volume_derivatives) {
        derivative *= total_feed_amount_mol;
    }
    result.local_variables = solution;
    result.local_variables.insert(
        result.local_variables.end(), multipliers.begin(), multipliers.end()
    );
    result.local_lower_bounds = lower;
    result.local_lower_bounds.insert(
        result.local_lower_bounds.end(), constraint_count, -1.0e19
    );
    result.local_upper_bounds = upper;
    result.local_upper_bounds.insert(
        result.local_upper_bounds.end(), constraint_count, 1.0e19
    );
    if (result.sensitivity_status == "available") {
        populate_outputs(
            result,
            evaluated,
            parameters,
            system.g_ref_parameter_derivatives,
            temperature_k,
            pressure_pa,
            total_feed_amount_mol,
            layouts,
            outputs
        );
    }
    if (result.sensitivity_status != "available") {
        result.accepted = false;
        result.failure_reason = "reactive_finite_phase_sensitivity_failed";
    }
    result.provider_evaluation_count = static_cast<std::size_t>(provider.evaluation_count());
    return result;
}

ReactiveFinitePhaseResult evaluate_provider_reactive_finite_phases_local(
    const CompiledReactionSystem& system,
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    const std::vector<std::array<double, 2>>& packing_fraction_intervals,
    const std::vector<double>& local_variables,
    const std::vector<double>& local_lower_bounds,
    const std::vector<double>& local_upper_bounds,
    const ProviderActiveParameterSet* active_parameters,
    const PhaseDirectionCompilation* phase_directions,
    const std::vector<std::vector<std::size_t>>& phase_support_indices,
    const std::vector<ReactiveFinitePhaseOutput>& outputs
) {
    const ProviderActiveParameterSet parameters = active_parameters == nullptr
        ? ProviderActiveParameterSet{}
        : *active_parameters;
    const std::vector<PhaseLayout> layouts = make_phase_layouts(
        system.charges, phase_support_indices
    );
    const std::size_t primal_count =
        layouts.back().offset + layouts.back().dimension();
    NormalizedBalances balances = normalized_balances(
        system, phase_directions, phase_support_indices
    );
    const double amount_scale = std::accumulate(
        system.feed_amounts.begin(), system.feed_amounts.end(), 0.0
    );
    for (double& coefficient : balances.matrix.values) {
        coefficient *= amount_scale;
    }
    const std::size_t constraint_count = balances.matrix.rows;
    if (local_variables.size() != primal_count + constraint_count
        || local_lower_bounds.size() != local_variables.size()
        || local_upper_bounds.size() != local_variables.size()) {
        throw std::invalid_argument("reactive finite-phase local dimensions changed");
    }
    const std::vector<double> primal(
        local_variables.begin(),
        local_variables.begin() + static_cast<std::ptrdiff_t>(primal_count)
    );
    const std::vector<double> multipliers(
        local_variables.begin() + static_cast<std::ptrdiff_t>(primal_count),
        local_variables.end()
    );
    const MultiphaseEvaluation evaluated = evaluate_multiphase(
        provider,
        parameters,
        system,
        balances,
        temperature_k,
        pressure_pa,
        packing_fraction_intervals,
        primal,
        multipliers,
        layouts
    );
    ReactiveFinitePhaseResult result;
    result.local_variables = local_variables;
    result.local_lower_bounds = local_lower_bounds;
    result.local_upper_bounds = local_upper_bounds;
    result.objective_over_rt = amount_scale * evaluated.objective;
    result.balance_inf_norm = inf_norm(evaluated.constraints);
    for (const PhaseBlock& block : evaluated.phases) {
        std::vector<double> amounts = block.amounts;
        for (double& amount : amounts) {
            amount *= amount_scale;
        }
        result.phase_amounts_mol.push_back(std::move(amounts));
        result.phase_volumes_m3.push_back(amount_scale * block.volume);
        result.phase_pressures_pa.push_back(block.pressure);
        result.phase_packing_fractions.push_back(block.packing);
        result.phase_chemical_potentials_over_rt.push_back(block.chemical_potentials);
        std::vector<double> log_fugacities = block.chemical_potentials;
        for (std::size_t species = 0; species < system.species_count; ++species) {
            log_fugacities[species] += std::log(
                kGasConstantJPerMolK * temperature_k
            ) - system.g_ref[species];
        }
        result.phase_log_fugacities.push_back(std::move(log_fugacities));
    }
    populate_sensitivities(
        result,
        evaluated,
        parameters,
        balances,
        system.g_ref_pressure_derivatives_per_pa,
        system.g_ref_parameter_derivatives,
        multipliers,
        temperature_k,
        system.species_count,
        layouts
    );
    for (double& derivative : result.amount_derivatives) {
        derivative *= amount_scale;
    }
    for (double& derivative : result.volume_derivatives) {
        derivative *= amount_scale;
    }
    for (double& derivative : result.input_amount_derivatives) {
        derivative *= amount_scale;
    }
    for (double& derivative : result.input_volume_derivatives) {
        derivative *= amount_scale;
    }
    populate_outputs(
        result,
        evaluated,
        parameters,
        system.g_ref_parameter_derivatives,
        temperature_k,
        pressure_pa,
        amount_scale,
        layouts,
        outputs
    );
    result.accepted = true;
    result.solver_status = "local_evaluation";
    result.numerical_status = "passed";
    result.physical_status = "not_adjudicated";
    result.provider_domain_status = "passed";
    result.provider_evaluation_count = static_cast<std::size_t>(provider.evaluation_count());
    return result;
}

}  // namespace epcsaft_engine::equilibrium
