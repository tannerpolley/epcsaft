#include "phase_equilibrium.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

#include <Eigen/Dense>
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpSolveStatistics.hpp>
#include <coin/IpTNLP.hpp>

namespace epcsaft_engine::equilibrium {
namespace {

constexpr double kGasConstantJPerMolK = 8.31446261815324;
constexpr double kPackingFractionMinimum = 1.0e-7;
constexpr double kPackingFractionMaximum = 0.70;
constexpr double kAmountFloorRelative = 1.0e-9;
constexpr double kAmountFloorAbsoluteMol = 1.0e-12;
constexpr double kMaterialToleranceMol = 2.0e-9;
constexpr double kChargeToleranceMol = 2.0e-9;
constexpr double kPressureRelativeTolerance = 1.0e-8;
constexpr double kChemicalPotentialTolerance = 2.0e-8;
constexpr double kStationarityTolerance = 2.0e-8;
constexpr double kPhaseDistinctTolerance = 2.0e-7;
constexpr double kRolePackingSeparation = 1.0e-5;
constexpr double kLiquidCompressibilityFactorMaximum = 0.25;
constexpr double kVaporCompressibilityFactorMinimum = 0.45;
constexpr double kReducedHessianTolerance = 2.0e-7;
constexpr double kKktConditionMaximum = 1.0e12;
constexpr int kStartBudgetPerTopology = 12;

using Matrix = Eigen::Matrix<
    double,
    Eigen::Dynamic,
    Eigen::Dynamic,
    Eigen::RowMajor
>;

bool finite(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
    });
}

std::vector<std::string> component_ids(const ProviderContext& provider) {
    const auto& transport = provider.transport();
    if (transport.component_count < 2 || transport.component_ids == nullptr) {
        throw std::invalid_argument("phase-equilibrium provider topology is incomplete");
    }
    std::vector<std::string> result;
    result.reserve(transport.component_count);
    for (std::size_t component = 0; component < transport.component_count; ++component) {
        const char* identity = transport.component_ids[component];
        if (identity == nullptr || identity[0] == '\0') {
            throw std::invalid_argument("phase-equilibrium component identity is incomplete");
        }
        result.emplace_back(identity);
    }
    return result;
}

std::vector<double> component_charges(const ProviderContext& provider) {
    const auto& transport = provider.transport();
    std::vector<double> result(transport.component_count, 0.0);
    if (transport.component_charges != nullptr) {
        for (std::size_t component = 0; component < transport.component_count; ++component) {
            result[component] = static_cast<double>(transport.component_charges[component]);
        }
    }
    return result;
}

bool charged_topology(const std::vector<double>& charges) {
    return std::any_of(charges.begin(), charges.end(), [](double value) {
        return value != 0.0;
    });
}

std::size_t charge_pivot(const std::vector<double>& charges) {
    const auto pivot = std::find_if(charges.begin(), charges.end(), [](double value) {
        return value != 0.0;
    });
    return pivot == charges.end()
        ? charges.size()
        : static_cast<std::size_t>(pivot - charges.begin());
}

std::size_t variable_index(
    std::size_t phase,
    std::size_t coordinate,
    std::size_t phase_dimension
) {
    return phase * phase_dimension + coordinate;
}

std::vector<std::size_t> material_components(
    std::size_t component_count,
    const std::vector<double>& charges
) {
    std::vector<std::size_t> result;
    result.reserve(component_count);
    const std::size_t pivot = charge_pivot(charges);
    if (pivot == component_count) {
        for (std::size_t component = 0; component < component_count; ++component) {
            result.push_back(component);
        }
        return result;
    }
    for (std::size_t component = 0; component < component_count; ++component) {
        if (component != pivot) {
            result.push_back(component);
        }
    }
    return result;
}

struct PhaseChart {
    std::vector<double> amounts;
    double log_packing_fraction = 0.0;
    std::vector<double> jacobian;
    std::vector<double> hessians;
};

std::vector<std::size_t> components_with_charge_sign(
    const std::vector<double>& charges,
    int sign
) {
    std::vector<std::size_t> result;
    for (std::size_t component = 0; component < charges.size(); ++component) {
        if ((sign > 0 && charges[component] > 0.0)
            || (sign < 0 && charges[component] < 0.0)) {
            result.push_back(component);
        }
    }
    return result;
}

std::vector<std::size_t> neutral_components(
    const std::vector<double>& charges
) {
    std::vector<std::size_t> result;
    for (std::size_t component = 0; component < charges.size(); ++component) {
        if (charges[component] == 0.0) {
            result.push_back(component);
        }
    }
    return result;
}

std::size_t chart_dimension(const std::vector<double>& charges) {
    return charged_topology(charges) ? charges.size() : charges.size() + 1;
}

std::vector<double> softmax_with_reference(
    const std::vector<double>& logits
) {
    const double maximum = std::max(
        0.0,
        logits.empty() ? 0.0 : *std::max_element(logits.begin(), logits.end())
    );
    std::vector<double> weights;
    weights.reserve(logits.size() + 1);
    double total = std::exp(-maximum);
    for (double logit : logits) {
        weights.push_back(std::exp(logit - maximum));
        total += weights.back();
    }
    weights.push_back(std::exp(-maximum));
    for (double& weight : weights) {
        weight /= total;
    }
    return weights;
}

void set_charge_group_chart(
    PhaseChart& chart,
    const std::vector<double>& charges,
    const std::vector<std::size_t>& components,
    double charge_equivalents,
    const std::vector<double>& logits,
    std::size_t charge_coordinate,
    std::size_t logit_offset,
    std::size_t phase_dimension
) {
    const auto weights = softmax_with_reference(logits);
    for (std::size_t item = 0; item < components.size(); ++item) {
        const std::size_t component = components[item];
        const double charge_magnitude = std::abs(charges[component]);
        chart.amounts[component] = charge_equivalents
            * weights[item] / charge_magnitude;
        chart.jacobian[component * phase_dimension + charge_coordinate] =
            weights[item] / charge_magnitude;
        for (std::size_t logit = 0; logit < logits.size(); ++logit) {
            const double delta_item = item == logit ? 1.0 : 0.0;
            const double first = weights[item]
                * (delta_item - weights[logit]);
            const std::size_t coordinate = logit_offset + logit;
            chart.jacobian[component * phase_dimension + coordinate] =
                charge_equivalents * first / charge_magnitude;
            chart.hessians[
                (component * phase_dimension + charge_coordinate)
                    * phase_dimension + coordinate
            ] = first / charge_magnitude;
            chart.hessians[
                (component * phase_dimension + coordinate)
                    * phase_dimension + charge_coordinate
            ] = first / charge_magnitude;
            for (std::size_t other = 0; other < logits.size(); ++other) {
                const double delta_item_other = item == other ? 1.0 : 0.0;
                const double delta_logit_other = logit == other ? 1.0 : 0.0;
                const double second = weights[item] * (
                    (delta_item - weights[logit])
                        * (delta_item_other - weights[other])
                    - weights[logit]
                        * (delta_logit_other - weights[other])
                );
                chart.hessians[
                    (component * phase_dimension + coordinate)
                        * phase_dimension + logit_offset + other
                ] = charge_equivalents * second / charge_magnitude;
            }
        }
    }
}

PhaseChart evaluate_phase_chart(
    const std::vector<double>& variables,
    const std::vector<double>& charges
) {
    const std::size_t component_count = charges.size();
    const std::size_t phase_dimension = chart_dimension(charges);
    if (variables.size() != phase_dimension) {
        throw std::invalid_argument("phase-equilibrium chart dimension changed");
    }
    PhaseChart result;
    result.amounts.assign(component_count, 0.0);
    result.jacobian.assign((component_count + 1) * phase_dimension, 0.0);
    result.hessians.assign(
        (component_count + 1) * phase_dimension * phase_dimension, 0.0
    );
    if (!charged_topology(charges)) {
        for (std::size_t component = 0; component < component_count; ++component) {
            result.amounts[component] = variables[component];
            result.jacobian[component * phase_dimension + component] = 1.0;
        }
        result.log_packing_fraction = variables.back();
        result.jacobian[component_count * phase_dimension + phase_dimension - 1] = 1.0;
        return result;
    }

    const auto neutral = neutral_components(charges);
    const auto positive = components_with_charge_sign(charges, 1);
    const auto negative = components_with_charge_sign(charges, -1);
    if (positive.empty() || negative.empty()) {
        throw std::invalid_argument(
            "charged phase-equilibrium topology needs both charge signs"
        );
    }
    std::size_t coordinate = 0;
    for (std::size_t component : neutral) {
        result.amounts[component] = variables[coordinate];
        result.jacobian[component * phase_dimension + coordinate] = 1.0;
        ++coordinate;
    }
    const std::size_t charge_coordinate = coordinate++;
    const double charge_equivalents = variables[charge_coordinate];
    const std::size_t positive_offset = coordinate;
    std::vector<double> positive_logits(
        variables.begin() + static_cast<std::ptrdiff_t>(coordinate),
        variables.begin() + static_cast<std::ptrdiff_t>(
            coordinate + positive.size() - 1
        )
    );
    coordinate += positive.size() - 1;
    const std::size_t negative_offset = coordinate;
    std::vector<double> negative_logits(
        variables.begin() + static_cast<std::ptrdiff_t>(coordinate),
        variables.begin() + static_cast<std::ptrdiff_t>(
            coordinate + negative.size() - 1
        )
    );
    coordinate += negative.size() - 1;
    if (coordinate + 1 != phase_dimension) {
        throw std::invalid_argument("charged phase-equilibrium chart is incomplete");
    }
    set_charge_group_chart(
        result,
        charges,
        positive,
        charge_equivalents,
        positive_logits,
        charge_coordinate,
        positive_offset,
        phase_dimension
    );
    set_charge_group_chart(
        result,
        charges,
        negative,
        charge_equivalents,
        negative_logits,
        charge_coordinate,
        negative_offset,
        phase_dimension
    );
    result.log_packing_fraction = variables.back();
    result.jacobian[component_count * phase_dimension + phase_dimension - 1] = 1.0;
    return result;
}

std::string topology_name(const std::vector<PhaseRole>& roles) {
    std::string result;
    for (PhaseRole role : roles) {
        result += role == PhaseRole::Liquid ? 'L' : 'V';
    }
    return result;
}

struct PhaseObjective {
    double value = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    std::vector<double> state_parameter_hessian;
    std::vector<double> volume_parameter_gradient;
    std::vector<double> volume_state_gradient;
    std::vector<double> amounts;
    double volume = 0.0;
    double packing_fraction = 0.0;
    MixturePhaseEvaluation provider;
};

MixturePhaseEvaluation evaluate_physical_phase(
    const ProviderContext& provider,
    bool charged,
    double temperature_k,
    const std::vector<double>& amounts,
    double volume
) {
    if (charged) {
        const auto& transport = provider.transport();
        if (transport.component_count != amounts.size()
            || transport.component_charges == nullptr) {
            throw std::invalid_argument(
                "charged phase-equilibrium provider topology is incomplete"
            );
        }
        const double total_amount = std::accumulate(
            amounts.begin(), amounts.end(), 0.0
        );
        double ionic_amount = 0.0;
        for (std::size_t component = 0; component < amounts.size(); ++component) {
            if (transport.component_charges[component] != 0) {
                ionic_amount += amounts[component];
            }
        }
        if (!std::isfinite(total_amount) || total_amount <= 0.0
            || !std::isfinite(ionic_amount)) {
            throw std::domain_error(
                "charged phase-equilibrium state has invalid amounts"
            );
        }
        if (std::isfinite(transport.total_ion_mole_fraction_max)
            && ionic_amount / total_amount
                > transport.total_ion_mole_fraction_max + 1.0e-12) {
            throw std::domain_error(
                "phase ion fraction is outside the model domain"
            );
        }
    }
    return charged
        ? provider.evaluate_electrolyte(temperature_k, amounts, volume)
        : provider.evaluate_mixture(temperature_k, amounts, volume);
}

PhaseObjective evaluate_phase_objective(
    const ProviderContext& provider,
    bool charged,
    double temperature_k,
    double pressure_over_rt,
    const ProviderActiveParameterSet& active,
    const std::vector<double>& amounts,
    double log_packing_fraction
) {
    const std::size_t component_count = amounts.size();
    const std::size_t state_count = component_count + 1;
    const std::size_t parameter_count = active.identities.size();
    const auto geometry = provider.evaluate_inverse_packing_geometry(
        temperature_k, amounts, log_packing_fraction, active
    );
    MixturePhaseEvaluation physical;
    std::vector<double> state_parameter;
    if (parameter_count == 0) {
        physical = evaluate_physical_phase(
            provider, charged, temperature_k, amounts, geometry.volume_m3
        );
    } else {
        auto parameterized = provider.evaluate_parameterized_phase(
            temperature_k,
            amounts,
            geometry.volume_m3,
            kPackingFractionMinimum,
            kPackingFractionMaximum,
            active
        );
        physical = std::move(parameterized.phase);
        state_parameter = std::move(parameterized.state_parameter_derivatives);
    }
    if (physical.gradient.size() != state_count
        || physical.hessian.size() != state_count * state_count
        || geometry.gradient.size() != state_count + parameter_count
        || geometry.hessian.size()
            != (state_count + parameter_count) * (state_count + parameter_count)
        || (parameter_count != 0
            && state_parameter.size() != state_count * parameter_count)) {
        throw std::invalid_argument("phase-equilibrium provider tensor dimensions changed");
    }

    const double volume_gradient = physical.gradient[component_count] + pressure_over_rt;
    PhaseObjective result;
    result.value = physical.value + pressure_over_rt * geometry.volume_m3;
    result.gradient.assign(state_count, 0.0);
    result.hessian.assign(state_count * state_count, 0.0);
    result.state_parameter_hessian.assign(state_count * parameter_count, 0.0);
    result.volume_parameter_gradient.assign(parameter_count, 0.0);
    result.volume_state_gradient.assign(
        geometry.gradient.begin(), geometry.gradient.begin() + state_count
    );

    for (std::size_t row = 0; row < state_count; ++row) {
        const double direct_gradient = row < component_count
            ? physical.gradient[row]
            : 0.0;
        result.gradient[row] = direct_gradient
            + volume_gradient * geometry.gradient[row];
        for (std::size_t column = 0; column < state_count; ++column) {
            double value = physical.hessian[
                component_count * state_count + component_count
            ] * geometry.gradient[row] * geometry.gradient[column]
                + volume_gradient * geometry.hessian[
                    row * (state_count + parameter_count) + column
                ];
            if (row < component_count) {
                value += physical.hessian[
                    row * state_count + component_count
                ] * geometry.gradient[column];
            }
            if (column < component_count) {
                value += physical.hessian[
                    component_count * state_count + column
                ] * geometry.gradient[row];
            }
            if (row < component_count && column < component_count) {
                value += physical.hessian[row * state_count + column];
            }
            result.hessian[row * state_count + column] = value;
        }
    }

    for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
        const std::size_t extended = state_count + parameter;
        const double volume_parameter = geometry.gradient[extended];
        result.volume_parameter_gradient[parameter] = volume_parameter;
        const double volume_parameter_coupling =
            physical.hessian[component_count * state_count + component_count]
                * volume_parameter
            + state_parameter[component_count * parameter_count + parameter];
        for (std::size_t row = 0; row < state_count; ++row) {
            double value = geometry.gradient[row] * volume_parameter_coupling
                + volume_gradient * geometry.hessian[
                    row * (state_count + parameter_count) + extended
                ];
            if (row < component_count) {
                value += state_parameter[row * parameter_count + parameter]
                    + physical.hessian[row * state_count + component_count]
                        * volume_parameter;
            }
            result.state_parameter_hessian[row * parameter_count + parameter] = value;
        }
    }
    result.amounts = amounts;
    result.volume = geometry.volume_m3;
    result.packing_fraction = std::exp(log_packing_fraction);
    result.provider = std::move(physical);
    return result;
}

PhaseObjective project_phase_objective(
    PhaseObjective full,
    const PhaseChart& chart,
    std::size_t parameter_count
) {
    const std::size_t full_dimension = chart.amounts.size() + 1;
    const std::size_t projected_dimension = chart.jacobian.size() / full_dimension;
    PhaseObjective result = full;
    result.gradient.assign(projected_dimension, 0.0);
    result.hessian.assign(projected_dimension * projected_dimension, 0.0);
    result.state_parameter_hessian.assign(
        projected_dimension * parameter_count, 0.0
    );
    result.volume_state_gradient.assign(projected_dimension, 0.0);
    for (std::size_t row = 0; row < projected_dimension; ++row) {
        for (std::size_t state = 0; state < full_dimension; ++state) {
            const double derivative = chart.jacobian[
                state * projected_dimension + row
            ];
            result.gradient[row] += derivative * full.gradient[state];
            result.volume_state_gradient[row] +=
                derivative * full.volume_state_gradient[state];
            for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
                result.state_parameter_hessian[row * parameter_count + parameter] +=
                    derivative * full.state_parameter_hessian[
                        state * parameter_count + parameter
                    ];
            }
        }
        for (std::size_t column = 0; column < projected_dimension; ++column) {
            double value = 0.0;
            for (std::size_t left = 0; left < full_dimension; ++left) {
                for (std::size_t right = 0; right < full_dimension; ++right) {
                    value += chart.jacobian[left * projected_dimension + row]
                        * full.hessian[left * full_dimension + right]
                        * chart.jacobian[right * projected_dimension + column];
                }
                value += full.gradient[left] * chart.hessians[
                    (left * projected_dimension + row) * projected_dimension
                        + column
                ];
            }
            result.hessian[row * projected_dimension + column] = value;
        }
    }
    return result;
}

struct ProblemEvaluation {
    double objective = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    std::vector<double> constraints;
    std::vector<double> jacobian;
    std::vector<double> constraint_hessians;
    std::vector<PhaseObjective> phases;
};

struct ConstraintEvaluation {
    std::vector<double> values;
    std::vector<double> jacobian;
    std::vector<double> hessians;
};

ConstraintEvaluation evaluate_constraints(
    const PhaseEquilibriumProblem& problem,
    const std::vector<PhaseRole>& roles,
    const std::vector<double>& charges,
    const std::vector<double>& variables
) {
    const std::size_t component_count = problem.feed_amounts_mol.size();
    const std::size_t phase_count = roles.size();
    const std::size_t phase_dimension = chart_dimension(charges);
    const std::size_t dimension = phase_count * phase_dimension;
    const auto balanced_components = material_components(component_count, charges);
    if (variables.size() != dimension) {
        throw std::invalid_argument("phase-equilibrium variable dimension changed");
    }
    ConstraintEvaluation result;
    result.values.assign(balanced_components.size(), 0.0);
    result.jacobian.assign(balanced_components.size() * dimension, 0.0);
    result.hessians.assign(
        balanced_components.size() * dimension * dimension, 0.0
    );
    std::vector<PhaseChart> charts;
    charts.reserve(phase_count);
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        const std::size_t offset = phase * phase_dimension;
        charts.push_back(evaluate_phase_chart(
            std::vector<double>(
                variables.begin() + static_cast<std::ptrdiff_t>(offset),
                variables.begin()
                    + static_cast<std::ptrdiff_t>(offset + phase_dimension)
            ),
            charges
        ));
    }
    for (std::size_t row = 0; row < balanced_components.size(); ++row) {
        const std::size_t component = balanced_components[row];
        result.values[row] = -problem.feed_amounts_mol[component];
        for (std::size_t phase = 0; phase < phase_count; ++phase) {
            const std::size_t offset = phase * phase_dimension;
            result.values[row] += charts[phase].amounts[component];
            for (std::size_t column = 0; column < phase_dimension; ++column) {
                result.jacobian[row * dimension + offset + column] =
                    charts[phase].jacobian[
                        component * phase_dimension + column
                    ];
                for (std::size_t other = 0; other < phase_dimension; ++other) {
                    result.hessians[
                        (row * dimension + offset + column) * dimension
                            + offset + other
                    ] = charts[phase].hessians[
                        (component * phase_dimension + column) * phase_dimension
                            + other
                    ];
                }
            }
        }
    }
    return result;
}

ProblemEvaluation evaluate_problem(
    const ProviderContext& provider,
    const PhaseEquilibriumProblem& problem,
    const std::vector<PhaseRole>& roles,
    const std::vector<double>& charges,
    const std::vector<double>& variables
) {
    const std::size_t phase_count = roles.size();
    const std::size_t phase_dimension = chart_dimension(charges);
    const std::size_t dimension = phase_count * phase_dimension;
    const bool charged = charged_topology(charges);
    if (variables.size() != dimension) {
        throw std::invalid_argument("phase-equilibrium variable dimension changed");
    }
    ProblemEvaluation result;
    result.gradient.assign(dimension, 0.0);
    result.hessian.assign(dimension * dimension, 0.0);
    ConstraintEvaluation constraints = evaluate_constraints(
        problem, roles, charges, variables
    );
    result.constraints = std::move(constraints.values);
    result.jacobian = std::move(constraints.jacobian);
    result.constraint_hessians = std::move(constraints.hessians);
    result.phases.reserve(phase_count);
    const double pressure_over_rt = problem.pressure_pa
        / (kGasConstantJPerMolK * problem.temperature_k);
    const ProviderActiveParameterSet active = provider.active_parameters();
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        const std::size_t offset = phase * phase_dimension;
        std::vector<double> chart_variables(
            variables.begin() + static_cast<std::ptrdiff_t>(offset),
            variables.begin() + static_cast<std::ptrdiff_t>(offset + phase_dimension)
        );
        PhaseChart chart = evaluate_phase_chart(chart_variables, charges);
        PhaseObjective block = project_phase_objective(evaluate_phase_objective(
            provider,
            charged,
            problem.temperature_k,
            pressure_over_rt,
            active,
            chart.amounts,
            chart.log_packing_fraction
        ), chart, active.identities.size());
        result.objective += block.value;
        for (std::size_t row = 0; row < phase_dimension; ++row) {
            result.gradient[offset + row] = block.gradient[row];
            for (std::size_t column = 0; column < phase_dimension; ++column) {
                result.hessian[(offset + row) * dimension + offset + column] =
                    block.hessian[row * phase_dimension + column];
            }
        }
        result.phases.push_back(std::move(block));
    }
    return result;
}

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
        case Ipopt::Invalid_Number_Detected:
            return "invalid_number_detected";
        default:
            return "ipopt_status_" + std::to_string(static_cast<int>(status));
    }
}

struct Start {
    std::string identity;
    std::vector<double> variables;
};

struct Attempt {
    bool converged = false;
    std::string solver_status = "not_evaluated";
    std::string callback_error;
    int iterations = 0;
    double constraint_violation = std::numeric_limits<double>::infinity();
    std::vector<PhaseRole> roles;
    std::vector<double> variables;
    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;
    std::vector<double> multipliers;
    std::vector<double> lower_multipliers;
    std::vector<double> upper_multipliers;
};

class PhaseTnlp final : public Ipopt::TNLP {
public:
    PhaseTnlp(
        const ProviderContext& provider,
        const PhaseEquilibriumProblem& problem,
        std::vector<PhaseRole> roles,
        std::vector<double> charges,
        std::vector<double> initial
    )
        : provider_(provider),
          problem_(problem),
          roles_(std::move(roles)),
          charges_(std::move(charges)),
          initial_(std::move(initial)) {
        const std::size_t component_count = problem_.feed_amounts_mol.size();
        const std::size_t phase_dimension = chart_dimension(charges_);
        const std::size_t dimension = roles_.size() * phase_dimension;
        lower_.assign(dimension, 0.0);
        upper_.assign(dimension, 0.0);
        for (std::size_t phase = 0; phase < roles_.size(); ++phase) {
            std::size_t coordinate = 0;
            for (std::size_t component : neutral_components(charges_)) {
                const std::size_t index = variable_index(
                    phase, coordinate, phase_dimension
                );
                lower_[index] = std::max(
                    kAmountFloorAbsoluteMol,
                    kAmountFloorRelative * problem_.feed_amounts_mol[component]
                );
                upper_[index] = problem_.feed_amounts_mol[component] + 1.0;
                ++coordinate;
            }
            if (!charged_topology(charges_)) {
                for (; coordinate < component_count; ++coordinate) {
                    const std::size_t index = variable_index(
                        phase, coordinate, phase_dimension
                    );
                    lower_[index] = std::max(
                        kAmountFloorAbsoluteMol,
                        kAmountFloorRelative * problem_.feed_amounts_mol[coordinate]
                    );
                    upper_[index] = problem_.feed_amounts_mol[coordinate] + 1.0;
                }
            } else {
                double positive_charge_equivalents = 0.0;
                for (std::size_t component = 0; component < component_count; ++component) {
                    if (charges_[component] > 0.0) {
                        positive_charge_equivalents += charges_[component]
                            * problem_.feed_amounts_mol[component];
                    }
                }
                const std::size_t charge_index = variable_index(
                    phase, coordinate++, phase_dimension
                );
                lower_[charge_index] = std::max(
                    kAmountFloorAbsoluteMol,
                    kAmountFloorRelative * positive_charge_equivalents
                );
                upper_[charge_index] = positive_charge_equivalents + 1.0;
                const std::size_t logit_count =
                    components_with_charge_sign(charges_, 1).size()
                    + components_with_charge_sign(charges_, -1).size() - 2;
                for (std::size_t logit = 0; logit < logit_count; ++logit) {
                    const std::size_t index = variable_index(
                        phase, coordinate++, phase_dimension
                    );
                    lower_[index] = -30.0;
                    upper_[index] = 30.0;
                }
            }
            lower_[variable_index(phase, phase_dimension - 1, phase_dimension)] =
                std::log(kPackingFractionMinimum);
            upper_[variable_index(phase, phase_dimension - 1, phase_dimension)] =
                std::log(kPackingFractionMaximum);
        }
        solution_.assign(dimension, 0.0);
        lower_multipliers_.assign(dimension, 0.0);
        upper_multipliers_.assign(dimension, 0.0);
        const std::size_t constraints = material_components(
            component_count, charges_
        ).size();
        multipliers_.assign(constraints, 0.0);
    }

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
        if (n != static_cast<Ipopt::Index>(initial_.size())
            || m != static_cast<Ipopt::Index>(multipliers_.size())) {
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
        bool init_z,
        Ipopt::Number*,
        Ipopt::Number*,
        Ipopt::Index m,
        bool init_lambda,
        Ipopt::Number*
    ) override {
        if (n != static_cast<Ipopt::Index>(initial_.size())
            || m != static_cast<Ipopt::Index>(multipliers_.size())
            || !init_x || init_z || init_lambda) {
            return false;
        }
        std::copy(initial_.begin(), initial_.end(), x);
        return true;
    }

    bool eval_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number& objective
    ) override {
        try {
            objective = evaluate(n, x).objective;
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            return false;
        }
    }

    bool eval_grad_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number* gradient
    ) override {
        try {
            const auto value = evaluate(n, x);
            std::copy(value.gradient.begin(), value.gradient.end(), gradient);
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            return false;
        }
    }

    bool eval_g(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Index m,
        Ipopt::Number* constraints
    ) override {
        try {
            const auto value = evaluate_constraints(
                problem_, roles_, charges_, std::vector<double>(x, x + n)
            );
            if (m != static_cast<Ipopt::Index>(value.values.size())) {
                return false;
            }
            std::copy(value.values.begin(), value.values.end(), constraints);
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            return false;
        }
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
        if (nonzero_count != n * m) {
            return false;
        }
        if (values == nullptr) {
            for (Ipopt::Index row = 0; row < m; ++row) {
                for (Ipopt::Index column = 0; column < n; ++column) {
                    rows[row * n + column] = row;
                    columns[row * n + column] = column;
                }
            }
            return true;
        }
        try {
            const auto value = evaluate_constraints(
                problem_, roles_, charges_, std::vector<double>(x, x + n)
            );
            std::copy(value.jacobian.begin(), value.jacobian.end(), values);
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            return false;
        }
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
        if (m != static_cast<Ipopt::Index>(multipliers_.size())
            || nonzero_count != n * (n + 1) / 2) {
            return false;
        }
        if (values == nullptr) {
            std::size_t index = 0;
            for (Ipopt::Index row = 0; row < n; ++row) {
                for (Ipopt::Index column = 0; column <= row; ++column) {
                    rows[index] = row;
                    columns[index] = column;
                    ++index;
                }
            }
            return true;
        }
        try {
            const auto value = evaluate(n, x);
            std::size_t index = 0;
            for (Ipopt::Index row = 0; row < n; ++row) {
                for (Ipopt::Index column = 0; column <= row; ++column) {
                    double hessian = objective_factor
                        * value.hessian[static_cast<std::size_t>(row * n + column)];
                    for (Ipopt::Index constraint = 0; constraint < m; ++constraint) {
                        hessian += multipliers[constraint]
                            * value.constraint_hessians[
                                static_cast<std::size_t>(
                                    (constraint * n + row) * n + column
                                )
                            ];
                    }
                    values[index++] = hessian;
                }
            }
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            return false;
        }
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
            && m == static_cast<Ipopt::Index>(multipliers_.size())
            && x != nullptr && z_l != nullptr && z_u != nullptr
            && lambda != nullptr) {
            std::copy(x, x + n, solution_.begin());
            std::copy(z_l, z_l + n, lower_multipliers_.begin());
            std::copy(z_u, z_u + n, upper_multipliers_.begin());
            std::copy(lambda, lambda + m, multipliers_.begin());
        }
    }

    [[nodiscard]] const std::vector<double>& solution() const { return solution_; }
    [[nodiscard]] const std::vector<double>& lower_bounds() const { return lower_; }
    [[nodiscard]] const std::vector<double>& upper_bounds() const { return upper_; }
    [[nodiscard]] const std::vector<double>& multipliers() const { return multipliers_; }
    [[nodiscard]] const std::vector<double>& lower_multipliers() const {
        return lower_multipliers_;
    }
    [[nodiscard]] const std::vector<double>& upper_multipliers() const {
        return upper_multipliers_;
    }
    [[nodiscard]] const std::string& callback_error() const { return callback_error_; }

private:
    ProblemEvaluation evaluate(Ipopt::Index n, const Ipopt::Number* x) const {
        if (n != static_cast<Ipopt::Index>(initial_.size())) {
            throw std::invalid_argument("phase-equilibrium NLP dimension changed");
        }
        return evaluate_problem(
            provider_,
            problem_,
            roles_,
            charges_,
            std::vector<double>(x, x + n)
        );
    }

    const ProviderContext& provider_;
    PhaseEquilibriumProblem problem_;
    std::vector<PhaseRole> roles_;
    std::vector<double> charges_;
    std::vector<double> initial_;
    std::vector<double> lower_;
    std::vector<double> upper_;
    std::vector<double> solution_;
    std::vector<double> multipliers_;
    std::vector<double> lower_multipliers_;
    std::vector<double> upper_multipliers_;
    std::string callback_error_;
};

Attempt run_ipopt(
    const ProviderContext& provider,
    const PhaseEquilibriumProblem& problem,
    const std::vector<PhaseRole>& roles,
    const std::vector<double>& charges,
    const Start& start
) {
    auto* raw_problem = new PhaseTnlp(
        provider, problem, roles, charges, start.variables
    );
    Ipopt::SmartPtr<Ipopt::TNLP> nlp = raw_problem;
    Ipopt::SmartPtr<Ipopt::IpoptApplication> application = IpoptApplicationFactory();
    application->Options()->SetIntegerValue("print_level", 0);
    application->Options()->SetStringValue("sb", "yes");
    application->Options()->SetIntegerValue("max_iter", 400);
    application->Options()->SetNumericValue("tol", 1.0e-12);
    application->Options()->SetNumericValue("acceptable_tol", 1.0e-11);
    application->Options()->SetIntegerValue("acceptable_iter", 0);
    application->Options()->SetNumericValue("constr_viol_tol", 1.0e-12);
    application->Options()->SetStringValue("jacobian_approximation", "exact");
    application->Options()->SetStringValue("hessian_approximation", "exact");
    application->Options()->SetStringValue("nlp_scaling_method", "gradient-based");
    application->Options()->SetNumericValue("bound_relax_factor", 0.0);
    application->Options()->SetStringValue("honor_original_bounds", "yes");
    application->Options()->SetStringValue("check_derivatives_for_naninf", "yes");

    Attempt result;
    result.roles = roles;
    const auto initialize_status = application->Initialize();
    if (initialize_status != Ipopt::Solve_Succeeded) {
        result.solver_status = "initialization_" + ipopt_status_name(initialize_status);
        return result;
    }
    const auto status = application->OptimizeTNLP(nlp);
    result.solver_status = ipopt_status_name(status);
    result.converged = status == Ipopt::Solve_Succeeded;
    result.variables = raw_problem->solution();
    result.lower_bounds = raw_problem->lower_bounds();
    result.upper_bounds = raw_problem->upper_bounds();
    result.multipliers = raw_problem->multipliers();
    result.lower_multipliers = raw_problem->lower_multipliers();
    result.upper_multipliers = raw_problem->upper_multipliers();
    result.callback_error = raw_problem->callback_error();
    const Ipopt::SmartPtr<Ipopt::SolveStatistics> statistics = application->Statistics();
    if (IsValid(statistics)) {
        result.iterations = statistics->IterationCount();
        double dual = 0.0;
        double complementarity = 0.0;
        double kkt = 0.0;
        statistics->Infeasibilities(
            dual, result.constraint_violation, complementarity, kkt
        );
    }
    return result;
}

double role_packing_seed(PhaseRole role, std::size_t ordinal) {
    static constexpr double liquid[] = {0.45, 0.25, 0.60};
    static constexpr double vapor[] = {1.0e-4, 1.0e-3, 1.0e-2};
    return role == PhaseRole::Liquid
        ? liquid[ordinal % 3]
        : vapor[ordinal % 3];
}

std::vector<std::vector<double>> composition_directions(
    const std::vector<double>& feed,
    const std::vector<double>& charges
) {
    const std::size_t component_count = feed.size();
    std::vector<std::vector<double>> result;
    const std::size_t pivot = charge_pivot(charges);
    if (pivot == component_count) {
        for (std::size_t component = 0; component + 1 < component_count; ++component) {
            std::vector<double> direction(component_count, 0.0);
            direction[component] = 1.0;
            direction.back() = -1.0;
            result.push_back(std::move(direction));
        }
    } else {
        for (std::size_t component = 0; component < component_count; ++component) {
            if (component == pivot) {
                continue;
            }
            std::vector<double> direction(component_count, 0.0);
            direction[component] = 1.0;
            direction[pivot] = -charges[component] / charges[pivot];
            result.push_back(std::move(direction));
        }
    }
    return result;
}

std::vector<double> chart_variables_from_amounts(
    const std::vector<double>& amounts,
    const std::vector<double>& charges,
    double log_packing_fraction
) {
    std::vector<double> result;
    result.reserve(chart_dimension(charges));
    if (!charged_topology(charges)) {
        result = amounts;
        result.push_back(log_packing_fraction);
        return result;
    }
    for (std::size_t component : neutral_components(charges)) {
        result.push_back(amounts[component]);
    }
    const auto positive = components_with_charge_sign(charges, 1);
    const auto negative = components_with_charge_sign(charges, -1);
    double charge_equivalents = 0.0;
    for (std::size_t component : positive) {
        charge_equivalents += charges[component] * amounts[component];
    }
    if (!std::isfinite(charge_equivalents) || charge_equivalents <= 0.0) {
        throw std::invalid_argument("charged phase start has no positive charge");
    }
    result.push_back(charge_equivalents);
    for (const auto& group : {positive, negative}) {
        const std::size_t reference = group.back();
        const double reference_equivalents =
            std::abs(charges[reference]) * amounts[reference];
        if (reference_equivalents <= 0.0) {
            throw std::invalid_argument("charged phase start is not interior");
        }
        for (std::size_t item = 0; item + 1 < group.size(); ++item) {
            const std::size_t component = group[item];
            result.push_back(std::log(
                std::abs(charges[component]) * amounts[component]
                    / reference_equivalents
            ));
        }
    }
    result.push_back(log_packing_fraction);
    if (result.size() != chart_dimension(charges)) {
        throw std::invalid_argument("charged phase start chart is incomplete");
    }
    return result;
}

bool append_perturbed_start(
    std::vector<Start>& starts,
    const std::vector<double>& feed,
    const std::vector<double>& charges,
    const std::vector<PhaseRole>& roles,
    const std::vector<std::vector<double>>& directions,
    const std::vector<std::vector<double>>& coefficients,
    const std::string& identity
) {
    const std::size_t phase_count = roles.size();
    const std::size_t component_count = feed.size();
    std::vector<std::vector<double>> phase_amounts(
        phase_count,
        std::vector<double>(component_count, 0.0)
    );
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        for (std::size_t component = 0; component < component_count; ++component) {
            phase_amounts[phase][component] =
                feed[component] / static_cast<double>(phase_count);
        }
    }
    double amplitude = std::numeric_limits<double>::infinity();
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        for (std::size_t component = 0; component < component_count; ++component) {
            double combined = 0.0;
            for (std::size_t direction = 0; direction < directions.size(); ++direction) {
                combined += coefficients[phase][direction]
                    * directions[direction][component];
            }
            if (combined != 0.0) {
                amplitude = std::min(
                    amplitude,
                    0.35 * feed[component]
                        / (static_cast<double>(phase_count) * std::abs(combined))
                );
            }
        }
    }
    if (!std::isfinite(amplitude) || amplitude <= 0.0) {
        return false;
    }
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        for (std::size_t component = 0; component < component_count; ++component) {
            double combined = 0.0;
            for (std::size_t direction = 0; direction < directions.size(); ++direction) {
                combined += coefficients[phase][direction]
                    * directions[direction][component];
            }
            phase_amounts[phase][component] += amplitude * combined;
            if (phase_amounts[phase][component] <= 0.0) {
                return false;
            }
        }
    }
    std::vector<double> variables;
    variables.reserve(phase_count * chart_dimension(charges));
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        const auto block = chart_variables_from_amounts(
            phase_amounts[phase],
            charges,
            std::log(role_packing_seed(roles[phase], phase))
        );
        variables.insert(variables.end(), block.begin(), block.end());
    }
    starts.push_back({identity, std::move(variables)});
    return true;
}

std::vector<Start> make_starts(
    const std::vector<double>& feed,
    const std::vector<double>& charges,
    const std::vector<PhaseRole>& roles
) {
    const std::size_t phase_count = roles.size();
    const std::size_t component_count = feed.size();
    const std::size_t phase_dimension = chart_dimension(charges);
    std::vector<Start> starts;
    std::vector<double> base(phase_count * phase_dimension, 0.0);
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        std::vector<double> amounts(component_count, 0.0);
        for (std::size_t component = 0; component < component_count; ++component) {
            amounts[component] = feed[component] / static_cast<double>(phase_count);
        }
        const auto block = chart_variables_from_amounts(
            amounts,
            charges,
            std::log(role_packing_seed(roles[phase], phase))
        );
        std::copy(
            block.begin(),
            block.end(),
            base.begin() + static_cast<std::ptrdiff_t>(phase * phase_dimension)
        );
    }
    starts.push_back({"uniform", base});
    if (phase_count == 1) {
        for (std::size_t seed = 1; seed < 3; ++seed) {
            auto varied = base;
            varied[phase_dimension - 1] = std::log(role_packing_seed(roles[0], seed));
            starts.push_back({"density_seed_" + std::to_string(seed), std::move(varied)});
        }
        return starts;
    }

    const auto directions = composition_directions(feed, charges);
    for (std::size_t direction = 0;
         direction < directions.size() && starts.size() < kStartBudgetPerTopology;
         ++direction) {
        for (double sign : {1.0, -1.0}) {
            std::vector<std::vector<double>> coefficients(
                phase_count, std::vector<double>(directions.size(), 0.0)
            );
            coefficients[0][direction] = sign;
            coefficients[1][direction] = -sign;
            append_perturbed_start(
                starts,
                feed,
                charges,
                roles,
                directions,
                coefficients,
                "direction_" + std::to_string(direction)
                    + (sign > 0.0 ? "_positive" : "_negative")
            );
            if (starts.size() >= kStartBudgetPerTopology) {
                break;
            }
        }
    }
    if (phase_count == 3 && directions.size() >= 2
        && starts.size() < kStartBudgetPerTopology) {
        std::vector<std::vector<double>> coefficients(
            phase_count, std::vector<double>(directions.size(), 0.0)
        );
        coefficients[0][0] = 1.0;
        coefficients[1][1] = 1.0;
        coefficients[2][0] = -1.0;
        coefficients[2][1] = -1.0;
        append_perturbed_start(
            starts, feed, charges, roles, directions, coefficients, "three_phase_corner"
        );
        for (auto& phase : coefficients) {
            for (double& value : phase) {
                value = -value;
            }
        }
        append_perturbed_start(
            starts,
            feed,
            charges,
            roles,
            directions,
            coefficients,
            "three_phase_corner_reversed"
        );
    }
    if (starts.size() > kStartBudgetPerTopology) {
        starts.resize(kStartBudgetPerTopology);
    }
    return starts;
}

double hint_score(
    const Start& start,
    const std::vector<PhaseRole>& roles,
    const std::vector<PhaseHint>& hints,
    const std::vector<double>& hint_packing_fractions,
    const std::vector<double>& charges
) {
    if (hints.empty()) {
        return 0.0;
    }
    double score = 0.0;
    const std::size_t component_count = charges.size();
    const std::size_t phase_dimension = chart_dimension(charges);
    const std::size_t count = std::min(roles.size(), hints.size());
    for (std::size_t phase = 0; phase < count; ++phase) {
        if (roles[phase] != hints[phase].role) {
            score += 100.0;
        }
        const std::size_t offset = phase * phase_dimension;
        const PhaseChart chart = evaluate_phase_chart(
            std::vector<double>(
                start.variables.begin() + static_cast<std::ptrdiff_t>(offset),
                start.variables.begin()
                    + static_cast<std::ptrdiff_t>(offset + phase_dimension)
            ),
            charges
        );
        const double amount = std::accumulate(
            chart.amounts.begin(), chart.amounts.end(), 0.0
        );
        for (std::size_t component = 0; component < component_count; ++component) {
            score += std::abs(
                chart.amounts[component] / amount
                    - hints[phase].mole_fractions[component]
            );
        }
        score += 0.01 * std::abs(
            chart.log_packing_fraction
                - std::log(hint_packing_fractions[phase])
        );
    }
    return score;
}

int reorder_starts_from_hints(
    std::vector<Start>& starts,
    const ProviderContext& provider,
    const PhaseEquilibriumProblem& problem,
    const std::vector<PhaseRole>& roles,
    const std::vector<PhaseHint>& hints,
    const std::vector<double>& charges
) {
    if (hints.empty()) {
        return 0;
    }
    std::vector<double> hint_packing_fractions;
    hint_packing_fractions.reserve(hints.size());
    for (const auto& hint : hints) {
        hint_packing_fractions.push_back(
            provider.evaluate_packing_fraction(
                problem.temperature_k,
                hint.mole_fractions,
                1.0 / hint.molar_density_mol_m3
            ).value
        );
    }
    std::vector<std::string> original;
    original.reserve(starts.size());
    for (const auto& start : starts) {
        original.push_back(start.identity);
    }
    std::stable_sort(starts.begin(), starts.end(), [&](const Start& left, const Start& right) {
        return hint_score(
            left, roles, hints, hint_packing_fractions, charges
        ) < hint_score(
            right, roles, hints, hint_packing_fractions, charges
        );
    });
    int moved = 0;
    for (std::size_t index = 0; index < starts.size(); ++index) {
        moved += starts[index].identity != original[index] ? 1 : 0;
    }
    return moved;
}

struct Conditioning {
    std::size_t dimension = 0;
    std::size_t rank = 0;
    double condition = std::numeric_limits<double>::infinity();
    double minimum_reduced_eigenvalue = -std::numeric_limits<double>::infinity();
};

Matrix lagrangian_hessian(
    const ProblemEvaluation& evaluation,
    const std::vector<double>& multipliers
) {
    const std::size_t dimension = evaluation.gradient.size();
    if (multipliers.size() != evaluation.constraints.size()) {
        throw std::invalid_argument("phase-equilibrium multiplier dimension changed");
    }
    Eigen::Map<const Matrix> objective(
        evaluation.hessian.data(), dimension, dimension
    );
    Matrix result = objective;
    for (std::size_t constraint = 0; constraint < multipliers.size(); ++constraint) {
        Eigen::Map<const Matrix> block(
            evaluation.constraint_hessians.data()
                + static_cast<std::ptrdiff_t>(constraint * dimension * dimension),
            dimension,
            dimension
        );
        result += multipliers[constraint] * block;
    }
    return result;
}

Matrix equilibrated(Matrix matrix) {
    for (int pass = 0; pass < 4; ++pass) {
        for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
            const double norm = matrix.row(row).cwiseAbs().maxCoeff();
            if (norm > 0.0 && std::isfinite(norm)) {
                matrix.row(row) /= norm;
            }
        }
        for (Eigen::Index column = 0; column < matrix.cols(); ++column) {
            const double norm = matrix.col(column).cwiseAbs().maxCoeff();
            if (norm > 0.0 && std::isfinite(norm)) {
                matrix.col(column) /= norm;
            }
        }
    }
    return matrix;
}

Conditioning evaluate_conditioning(
    const ProblemEvaluation& evaluation,
    const std::vector<double>& multipliers
) {
    const std::size_t variable_count = evaluation.gradient.size();
    const std::size_t constraint_count = evaluation.constraints.size();
    Conditioning result;
    result.dimension = variable_count + constraint_count;
    const Matrix hessian = lagrangian_hessian(evaluation, multipliers);
    Eigen::Map<const Matrix> jacobian(
        evaluation.jacobian.data(), constraint_count, variable_count
    );
    Matrix kkt = Matrix::Zero(result.dimension, result.dimension);
    kkt.topLeftCorner(variable_count, variable_count) = hessian;
    kkt.topRightCorner(variable_count, constraint_count) = jacobian.transpose();
    kkt.bottomLeftCorner(constraint_count, variable_count) = jacobian;
    const Matrix scaled = equilibrated(kkt);
    Eigen::JacobiSVD<Matrix> decomposition(
        scaled, Eigen::ComputeFullU | Eigen::ComputeFullV
    );
    decomposition.setThreshold(
        4096.0 * std::numeric_limits<double>::epsilon()
            * static_cast<double>(result.dimension)
    );
    result.rank = static_cast<std::size_t>(decomposition.rank());
    const auto singular = decomposition.singularValues();
    if (result.rank == result.dimension && singular.size() != 0
        && singular[static_cast<Eigen::Index>(result.dimension - 1)] > 0.0) {
        result.condition = singular[0]
            / singular[static_cast<Eigen::Index>(result.dimension - 1)];
    }

    Eigen::JacobiSVD<Matrix> constraint_svd(
        jacobian, Eigen::ComputeFullV
    );
    constraint_svd.setThreshold(
        4096.0 * std::numeric_limits<double>::epsilon()
            * static_cast<double>(std::max(variable_count, constraint_count))
    );
    const Eigen::Index constraint_rank = constraint_svd.rank();
    const Eigen::Index nullity = static_cast<Eigen::Index>(variable_count)
        - constraint_rank;
    if (nullity == 0) {
        result.minimum_reduced_eigenvalue = 0.0;
    } else {
        const Matrix nullspace = constraint_svd.matrixV().rightCols(nullity);
        const Matrix reduced = nullspace.transpose() * hessian * nullspace;
        Eigen::SelfAdjointEigenSolver<Matrix> eigen(reduced);
        if (eigen.info() == Eigen::Success) {
            result.minimum_reduced_eigenvalue = eigen.eigenvalues().minCoeff();
        }
    }
    return result;
}

std::vector<PhaseEquilibriumState> states_from_evaluation(
    const ProblemEvaluation& evaluation,
    const std::vector<PhaseRole>& roles
) {
    std::vector<PhaseEquilibriumState> result;
    result.reserve(roles.size());
    for (std::size_t phase = 0; phase < roles.size(); ++phase) {
        const auto& block = evaluation.phases[phase];
        PhaseEquilibriumState state;
        state.role = roles[phase];
        state.amounts_mol = block.amounts;
        state.amount_mol = std::accumulate(
            state.amounts_mol.begin(), state.amounts_mol.end(), 0.0
        );
        state.mole_fractions.reserve(state.amounts_mol.size());
        for (double amount : state.amounts_mol) {
            state.mole_fractions.push_back(amount / state.amount_mol);
        }
        state.volume_m3 = block.volume;
        state.molar_density_mol_m3 = state.amount_mol / state.volume_m3;
        state.packing_fraction = block.packing_fraction;
        state.pressure_pa = block.provider.pressure_pa;
        state.chemical_potential_over_rt.assign(
            block.provider.gradient.begin(),
            block.provider.gradient.begin()
                + static_cast<std::ptrdiff_t>(state.amounts_mol.size())
        );
        result.push_back(std::move(state));
    }
    return result;
}

void canonicalize_attempt(
    Attempt& attempt,
    const PhaseEquilibriumProblem& problem,
    const std::vector<double>& charges,
    const ProviderContext& provider
) {
    if (attempt.variables.empty()) {
        return;
    }
    const std::size_t component_count = problem.feed_amounts_mol.size();
    const std::size_t phase_dimension = chart_dimension(charges);
    auto evaluation = evaluate_problem(
        provider, problem, attempt.roles, charges, attempt.variables
    );
    const auto ids = component_ids(provider);
    std::vector<std::size_t> identity_order(component_count);
    std::iota(identity_order.begin(), identity_order.end(), 0);
    std::stable_sort(
        identity_order.begin(), identity_order.end(),
        [&](std::size_t left, std::size_t right) { return ids[left] < ids[right]; }
    );
    std::vector<std::size_t> permutation(attempt.roles.size());
    std::iota(permutation.begin(), permutation.end(), 0);
    std::stable_sort(permutation.begin(), permutation.end(), [&](std::size_t left, std::size_t right) {
        if (attempt.roles[left] != attempt.roles[right]) {
            return attempt.roles[left] == PhaseRole::Liquid;
        }
        const double left_packing = evaluation.phases[left].packing_fraction;
        const double right_packing = evaluation.phases[right].packing_fraction;
        if (std::abs(left_packing - right_packing) > 1.0e-10) {
            return left_packing > right_packing;
        }
        const double left_amount = std::accumulate(
            evaluation.phases[left].amounts.begin(),
            evaluation.phases[left].amounts.end(),
            0.0
        );
        const double right_amount = std::accumulate(
            evaluation.phases[right].amounts.begin(),
            evaluation.phases[right].amounts.end(),
            0.0
        );
        if (std::abs(left_amount - right_amount) > 1.0e-10) {
            return left_amount > right_amount;
        }
        for (std::size_t component : identity_order) {
            const double left_fraction = evaluation.phases[left].amounts[component]
                / left_amount;
            const double right_fraction = evaluation.phases[right].amounts[component]
                / right_amount;
            if (std::abs(left_fraction - right_fraction) > 1.0e-10) {
                return left_fraction > right_fraction;
            }
        }
        return left < right;
    });
    bool identity = true;
    for (std::size_t index = 0; index < permutation.size(); ++index) {
        identity = identity && permutation[index] == index;
    }
    if (identity) {
        return;
    }
    auto reorder_blocks = [&](std::vector<double>& values) {
        std::vector<double> reordered(values.size(), 0.0);
        for (std::size_t target = 0; target < permutation.size(); ++target) {
            const std::size_t source = permutation[target];
            std::copy_n(
                values.begin() + static_cast<std::ptrdiff_t>(source * phase_dimension),
                phase_dimension,
                reordered.begin() + static_cast<std::ptrdiff_t>(target * phase_dimension)
            );
        }
        values = std::move(reordered);
    };
    reorder_blocks(attempt.variables);
    reorder_blocks(attempt.lower_bounds);
    reorder_blocks(attempt.upper_bounds);
    reorder_blocks(attempt.lower_multipliers);
    reorder_blocks(attempt.upper_multipliers);
    std::vector<PhaseRole> reordered_roles;
    reordered_roles.reserve(attempt.roles.size());
    for (std::size_t source : permutation) {
        reordered_roles.push_back(attempt.roles[source]);
    }
    attempt.roles = std::move(reordered_roles);
}

std::vector<std::string> phase_chart_identities(
    std::size_t phase,
    const std::vector<std::string>& ids,
    const std::vector<double>& charges
) {
    const std::string prefix = "phase[" + std::to_string(phase) + "].";
    std::vector<std::string> result;
    if (!charged_topology(charges)) {
        for (const auto& id : ids) {
            result.push_back(prefix + "amount[" + id + "]");
        }
        result.push_back(prefix + "log_packing_fraction");
        return result;
    }
    for (std::size_t component : neutral_components(charges)) {
        result.push_back(prefix + "amount[" + ids[component] + "]");
    }
    result.push_back(prefix + "charge_equivalents");
    for (const auto& group : {
             components_with_charge_sign(charges, 1),
             components_with_charge_sign(charges, -1),
         }) {
        const std::string& reference = ids[group.back()];
        for (std::size_t item = 0; item + 1 < group.size(); ++item) {
            result.push_back(
                prefix + "charge_log_ratio[" + ids[group[item]] + "/"
                    + reference + "]"
            );
        }
    }
    result.push_back(prefix + "log_packing_fraction");
    return result;
}

PhaseEquilibriumDiagnostics base_diagnostics() {
    PhaseEquilibriumDiagnostics result;
    result.material_balance_tolerance_mol = kMaterialToleranceMol;
    result.electroneutrality_tolerance_mol = kChargeToleranceMol;
    result.pressure_relative_tolerance = kPressureRelativeTolerance;
    result.chemical_potential_tolerance_over_rt = kChemicalPotentialTolerance;
    result.kkt_stationarity_tolerance = kStationarityTolerance;
    result.phase_distinct_tolerance = kPhaseDistinctTolerance;
    result.liquid_compressibility_factor_maximum =
        kLiquidCompressibilityFactorMaximum;
    result.vapor_compressibility_factor_minimum =
        kVaporCompressibilityFactorMinimum;
    result.reduced_hessian_relative_tolerance = kReducedHessianTolerance;
    result.kkt_condition_number_maximum = kKktConditionMaximum;
    return result;
}

PhaseEquilibriumDiagnostics certify_attempt(
    const PhaseEquilibriumProblem& problem,
    const std::vector<double>& charges,
    const Attempt& attempt,
    const ProblemEvaluation& evaluation,
    const std::vector<std::string>& ids
) {
    PhaseEquilibriumDiagnostics result = base_diagnostics();
    result.solver_status = attempt.solver_status;
    result.failure_stage = "local_equilibrium";
    const std::size_t component_count = problem.feed_amounts_mol.size();
    const std::size_t phase_count = attempt.roles.size();
    const std::size_t variable_count = evaluation.gradient.size();
    const auto balanced_components = material_components(component_count, charges);
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        for (const auto& identity : phase_chart_identities(phase, ids, charges)) {
            result.equation_state_identities.push_back(identity);
            result.equation_residual_identities.push_back(
                "stationarity." + identity
            );
        }
    }
    for (std::size_t component : balanced_components) {
        result.equation_state_identities.push_back(
            "multiplier.material_balance[" + ids[component] + "]"
        );
        result.equation_residual_identities.push_back(
            "material_balance[" + ids[component] + "]"
        );
    }

    std::vector<double> reconstructed(component_count, 0.0);
    for (const auto& phase : evaluation.phases) {
        for (std::size_t component = 0; component < component_count; ++component) {
            reconstructed[component] += phase.amounts[component];
        }
    }
    for (std::size_t component = 0; component < component_count; ++component) {
        result.material_balance_inf = std::max(
            result.material_balance_inf,
            std::abs(reconstructed[component] - problem.feed_amounts_mol[component])
        );
    }
    for (const auto& phase : evaluation.phases) {
        double charge = 0.0;
        for (std::size_t component = 0; component < component_count; ++component) {
            charge += charges[component] * phase.amounts[component];
        }
        result.electroneutrality_inf = std::max(
            result.electroneutrality_inf, std::abs(charge)
        );
        result.pressure_residual_inf = std::max(
            result.pressure_residual_inf,
            std::abs(phase.provider.pressure_pa / problem.pressure_pa - 1.0)
        );
    }
    if (phase_count > 1) {
        const std::size_t pivot = charge_pivot(charges);
        for (std::size_t phase = 1; phase < phase_count; ++phase) {
            for (std::size_t component : balanced_components) {
                double phase_potential =
                    evaluation.phases[phase].provider.gradient[component];
                double reference_potential =
                    evaluation.phases[0].provider.gradient[component];
                if (pivot != component_count) {
                    const double charge_ratio = charges[component] / charges[pivot];
                    phase_potential -= charge_ratio
                        * evaluation.phases[phase].provider.gradient[pivot];
                    reference_potential -= charge_ratio
                        * evaluation.phases[0].provider.gradient[pivot];
                }
                result.chemical_potential_residual_inf = std::max(
                    result.chemical_potential_residual_inf,
                    std::abs(phase_potential - reference_potential)
                );
            }
        }
    }
    for (std::size_t column = 0; column < variable_count; ++column) {
        double stationarity = evaluation.gradient[column]
            - attempt.lower_multipliers[column]
            + attempt.upper_multipliers[column];
        for (std::size_t row = 0; row < evaluation.constraints.size(); ++row) {
            stationarity += attempt.multipliers[row]
                * evaluation.jacobian[row * variable_count + column];
        }
        result.kkt_stationarity_inf = std::max(
            result.kkt_stationarity_inf, std::abs(stationarity)
        );
    }
    const Matrix certified_hessian = lagrangian_hessian(
        evaluation, attempt.multipliers
    );
    const Conditioning conditioning = evaluate_conditioning(
        evaluation, attempt.multipliers
    );
    result.minimum_reduced_hessian_eigenvalue =
        conditioning.minimum_reduced_eigenvalue;
    result.kkt_dimension = conditioning.dimension;
    result.kkt_rank = conditioning.rank;
    result.scaled_kkt_condition_number_2 = conditioning.condition;

    bool distinct = true;
    for (std::size_t left = 0; left < phase_count; ++left) {
        const double left_total = std::accumulate(
            evaluation.phases[left].amounts.begin(),
            evaluation.phases[left].amounts.end(),
            0.0
        );
        for (std::size_t right = left + 1; right < phase_count; ++right) {
            const double right_total = std::accumulate(
                evaluation.phases[right].amounts.begin(),
                evaluation.phases[right].amounts.end(),
                0.0
            );
            double composition_distance = 0.0;
            for (std::size_t component = 0; component < component_count; ++component) {
                composition_distance = std::max(
                    composition_distance,
                    std::abs(
                        evaluation.phases[left].amounts[component] / left_total
                            - evaluation.phases[right].amounts[component] / right_total
                    )
                );
            }
            const double packing_distance = std::abs(
                evaluation.phases[left].packing_fraction
                    - evaluation.phases[right].packing_fraction
            );
            distinct = distinct
                && (composition_distance > kPhaseDistinctTolerance
                    || packing_distance > kPhaseDistinctTolerance);
        }
    }
    bool role_separated = true;
    for (std::size_t liquid = 0; liquid < phase_count; ++liquid) {
        if (attempt.roles[liquid] != PhaseRole::Liquid) {
            continue;
        }
        for (std::size_t vapor = 0; vapor < phase_count; ++vapor) {
            if (attempt.roles[vapor] == PhaseRole::Vapor) {
                role_separated = role_separated
                    && evaluation.phases[liquid].packing_fraction
                        > evaluation.phases[vapor].packing_fraction
                            + kRolePackingSeparation;
            }
        }
    }
    bool roles_classified = true;
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        const double amount = std::accumulate(
            evaluation.phases[phase].amounts.begin(),
            evaluation.phases[phase].amounts.end(),
            0.0
        );
        const double compressibility_factor = problem.pressure_pa
            * evaluation.phases[phase].volume
            / (amount * kGasConstantJPerMolK * problem.temperature_k);
        roles_classified = roles_classified
            && (attempt.roles[phase] == PhaseRole::Liquid
                ? compressibility_factor <= kLiquidCompressibilityFactorMaximum
                : compressibility_factor >= kVaporCompressibilityFactorMinimum);
    }
    bool interior = attempt.lower_bounds.size() == variable_count
        && attempt.upper_bounds.size() == variable_count
        && attempt.lower_multipliers.size() == variable_count
        && attempt.upper_multipliers.size() == variable_count;
    for (std::size_t coordinate = 0; coordinate < variable_count && interior;
         ++coordinate) {
        const double scale = std::max({
            1.0,
            std::abs(attempt.variables[coordinate]),
            std::abs(attempt.lower_bounds[coordinate]),
            std::abs(attempt.upper_bounds[coordinate]),
        });
        const double margin = 4096.0
            * std::numeric_limits<double>::epsilon() * scale;
        interior = attempt.variables[coordinate]
                > attempt.lower_bounds[coordinate] + margin
            && attempt.variables[coordinate]
                < attempt.upper_bounds[coordinate] - margin
            && std::abs(attempt.lower_multipliers[coordinate])
                <= kStationarityTolerance
            && std::abs(attempt.upper_multipliers[coordinate])
                <= kStationarityTolerance;
    }
    const std::size_t phase_dimension = chart_dimension(charges);
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        for (std::size_t component = 0; component < component_count; ++component) {
            const double floor = std::max(
                kAmountFloorAbsoluteMol,
                kAmountFloorRelative * problem.feed_amounts_mol[component]
            );
            interior = interior
                && evaluation.phases[phase].amounts[component] > 2.0 * floor;
        }
        const double log_packing = attempt.variables[
            variable_index(phase, phase_dimension - 1, phase_dimension)
        ];
        interior = interior
            && log_packing > std::log(kPackingFractionMinimum) + 1.0e-8
            && log_packing < std::log(kPackingFractionMaximum) - 1.0e-8;
    }
    const double hessian_scale = std::max(
        1.0, certified_hessian.cwiseAbs().maxCoeff()
    );
    const bool accepted = attempt.converged && attempt.callback_error.empty()
        && std::isfinite(attempt.constraint_violation)
        && attempt.constraint_violation <= kMaterialToleranceMol
        && result.material_balance_inf <= kMaterialToleranceMol
        && result.electroneutrality_inf <= kChargeToleranceMol
        && result.pressure_residual_inf <= kPressureRelativeTolerance
        && result.chemical_potential_residual_inf <= kChemicalPotentialTolerance
        && result.kkt_stationarity_inf <= kStationarityTolerance
        && distinct && role_separated && roles_classified && interior
        && conditioning.rank == conditioning.dimension
        && std::isfinite(conditioning.condition)
        && conditioning.condition <= kKktConditionMaximum
        && std::isfinite(conditioning.minimum_reduced_eigenvalue)
        && conditioning.minimum_reduced_eigenvalue
            >= -kReducedHessianTolerance * hessian_scale;
    if (accepted) {
        result.numerical_status = "certified";
        result.physical_status = "certified";
        result.failure_stage.clear();
        return result;
    }
    result.numerical_status = attempt.converged ? "uncertified" : "not_converged";
    result.physical_status = "rejected";
    if (!attempt.callback_error.empty()) {
        result.failure_reason = "provider_evaluation_failed";
    } else if (!attempt.converged) {
        result.failure_reason = "local_solver_not_converged";
    } else if (!distinct) {
        result.failure_reason = "coalesced_phase_topology";
    } else if (!role_separated) {
        result.failure_reason = "declared_phase_roles_not_satisfied";
    } else if (!roles_classified) {
        result.failure_reason = "phase_role_classification_inconclusive";
    } else if (!interior) {
        result.failure_reason = "phase_boundary_or_domain_boundary";
    } else if (conditioning.rank != conditioning.dimension
        || !std::isfinite(conditioning.condition)
        || conditioning.condition > kKktConditionMaximum) {
        result.failure_reason = "singular_or_ill_conditioned_kkt_system";
    } else if (conditioning.minimum_reduced_eigenvalue
        < -kReducedHessianTolerance * hessian_scale) {
        result.failure_reason = "second_order_local_minimum_not_certified";
    } else {
        result.failure_reason = "equilibrium_residual_not_certified";
    }
    return result;
}

struct Candidate {
    Attempt attempt;
    ProblemEvaluation evaluation;
    PhaseEquilibriumDiagnostics diagnostics;
    std::vector<PhaseEquilibriumState> states;
    bool certified = false;
};

struct LinearSolve {
    std::size_t rank = 0;
    double condition = std::numeric_limits<double>::infinity();
    double backward_error = std::numeric_limits<double>::infinity();
    Matrix solution;
};

LinearSolve solve_implicit(const Matrix& matrix, const Matrix& partials) {
    LinearSolve result;
    Matrix scaled = matrix;
    Eigen::VectorXd row_scale = Eigen::VectorXd::Ones(matrix.rows());
    Eigen::VectorXd column_scale = Eigen::VectorXd::Ones(matrix.cols());
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
            * static_cast<double>(matrix.rows())
    );
    result.rank = static_cast<std::size_t>(decomposition.rank());
    const auto singular = decomposition.singularValues();
    if (result.rank == static_cast<std::size_t>(matrix.rows())
        && singular.size() != 0 && singular[singular.size() - 1] > 0.0) {
        result.condition = singular[0] / singular[singular.size() - 1];
    }
    if (result.rank != static_cast<std::size_t>(matrix.rows())
        || !std::isfinite(result.condition)
        || result.condition > kSolvedStateSensitivityConditionNumberMax) {
        return result;
    }
    const Matrix rhs = (-row_scale).asDiagonal() * partials;
    const Matrix scaled_solution = decomposition.solve(rhs);
    result.solution = column_scale.asDiagonal() * scaled_solution;
    const Matrix residual = matrix * result.solution + partials;
    const double numerator = residual.cwiseAbs().maxCoeff();
    const double denominator = std::max(
        1.0,
        matrix.cwiseAbs().maxCoeff() * result.solution.cwiseAbs().maxCoeff()
            + partials.cwiseAbs().maxCoeff()
    );
    result.backward_error = numerator / denominator;
    return result;
}

SolvedStateSensitivity evaluate_kkt_sensitivity(
    const ProviderContext& provider,
    const PhaseEquilibriumProblem& problem,
    const std::vector<double>& charges,
    const Candidate& candidate,
    const std::vector<std::string>& ids
) {
    SolvedStateSensitivity result;
    const ProviderActiveParameterSet active = provider.active_parameters();
    result.parameter_identities = active.identities;
    if (active.identities.empty()) {
        result.failure_reason = "no_active_parameters";
        return result;
    }
    const std::size_t component_count = problem.feed_amounts_mol.size();
    const std::size_t phase_count = candidate.attempt.roles.size();
    const std::size_t phase_dimension = chart_dimension(charges);
    const std::size_t variable_count = candidate.evaluation.gradient.size();
    const std::size_t constraint_count = candidate.evaluation.constraints.size();
    const std::size_t dimension = variable_count + constraint_count;
    const std::size_t parameter_count = active.identities.size();
    const auto balanced_components = material_components(component_count, charges);
    result.reference_phase_index = phase_count - 1;
    result.dimension = dimension;
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        for (const auto& identity : phase_chart_identities(phase, ids, charges)) {
            result.state_identities.push_back(identity);
            result.residual_identities.push_back("stationarity." + identity);
        }
        for (std::size_t component = 0; component < component_count; ++component) {
            result.output_identities.push_back(
                "phase[" + std::to_string(phase) + "].amount[" + ids[component] + "]"
            );
        }
        result.output_identities.push_back(
            "phase[" + std::to_string(phase) + "].volume_m3"
        );
    }
    for (std::size_t component : balanced_components) {
        result.state_identities.push_back(
            "multiplier.material_balance[" + ids[component] + "]"
        );
        result.residual_identities.push_back(
            "material_balance[" + ids[component] + "]"
        );
    }

    Matrix f_y = Matrix::Zero(dimension, dimension);
    const Matrix hessian = lagrangian_hessian(
        candidate.evaluation, candidate.attempt.multipliers
    );
    Eigen::Map<const Matrix> jacobian(
        candidate.evaluation.jacobian.data(), constraint_count, variable_count
    );
    f_y.topLeftCorner(variable_count, variable_count) = hessian;
    f_y.topRightCorner(variable_count, constraint_count) = jacobian.transpose();
    f_y.bottomLeftCorner(constraint_count, variable_count) = jacobian;
    Matrix f_p = Matrix::Zero(dimension, parameter_count);
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        const auto& block = candidate.evaluation.phases[phase];
        const std::size_t offset = phase * phase_dimension;
        for (std::size_t row = 0; row < phase_dimension; ++row) {
            for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
                f_p(offset + row, parameter) = block.state_parameter_hessian[
                    row * parameter_count + parameter
                ];
            }
        }
    }
    const LinearSolve solved = solve_implicit(f_y, f_p);
    result.rank = solved.rank;
    result.scaled_condition_number_2 = solved.condition;
    result.linear_solve_backward_error = solved.backward_error;
    result.residual_state_jacobian.assign(f_y.data(), f_y.data() + f_y.size());
    result.residual_parameter_jacobian.assign(f_p.data(), f_p.data() + f_p.size());
    if (solved.solution.rows() != static_cast<Eigen::Index>(dimension)
        || solved.solution.cols() != static_cast<Eigen::Index>(parameter_count)
        || solved.rank != dimension) {
        result.failure_reason = "singular_or_ill_conditioned_implicit_system";
        return result;
    }
    if (!std::isfinite(solved.backward_error)
        || solved.backward_error
            > solved_state_sensitivity_backward_error_max(dimension)
        || !solved.solution.allFinite()) {
        result.failure_reason =
            "implicit_linear_solve_backward_error_not_certified";
        return result;
    }
    result.state_parameter_jacobian.assign(
        solved.solution.data(), solved.solution.data() + solved.solution.size()
    );
    Matrix output = Matrix::Zero(
        phase_count * (component_count + 1), parameter_count
    );
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        const auto& block = candidate.evaluation.phases[phase];
        const std::size_t state_offset = phase * phase_dimension;
        const std::size_t output_offset = phase * (component_count + 1);
        const std::size_t chart_offset = phase * phase_dimension;
        const PhaseChart chart = evaluate_phase_chart(
            std::vector<double>(
                candidate.attempt.variables.begin()
                    + static_cast<std::ptrdiff_t>(chart_offset),
                candidate.attempt.variables.begin()
                    + static_cast<std::ptrdiff_t>(chart_offset + phase_dimension)
            ),
            charges
        );
        for (std::size_t component = 0; component < component_count; ++component) {
            for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
                for (std::size_t state = 0; state < phase_dimension; ++state) {
                    output(output_offset + component, parameter) += chart.jacobian[
                        component * phase_dimension + state
                    ] * solved.solution(state_offset + state, parameter);
                }
            }
        }
        for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
            double derivative = block.volume_parameter_gradient[parameter];
            for (std::size_t state = 0; state < phase_dimension; ++state) {
                derivative += block.volume_state_gradient[state]
                    * solved.solution(state_offset + state, parameter);
            }
            output(output_offset + component_count, parameter) = derivative;
        }
    }
    result.output_parameter_jacobian.assign(
        output.data(), output.data() + output.size()
    );
    result.status = "available";
    return result;
}

std::vector<std::vector<PhaseRole>> candidate_topologies(
    const PhaseEquilibriumProblem& problem
) {
    if (problem.mode == PhaseEquilibriumMode::Constrained) {
        return {problem.constrained_roles};
    }
    return {
        {PhaseRole::Liquid},
        {PhaseRole::Vapor},
        {PhaseRole::Liquid, PhaseRole::Liquid},
        {PhaseRole::Liquid, PhaseRole::Vapor},
        {PhaseRole::Liquid, PhaseRole::Liquid, PhaseRole::Liquid},
        {PhaseRole::Liquid, PhaseRole::Liquid, PhaseRole::Vapor},
    };
}

void validate_problem(
    const ProviderContext& provider,
    const PhaseEquilibriumProblem& problem,
    const std::vector<double>& charges
) {
    const auto& transport = provider.transport();
    const std::size_t component_count = transport.component_count;
    if (!std::isfinite(problem.temperature_k) || problem.temperature_k <= 0.0
        || !std::isfinite(problem.pressure_pa) || problem.pressure_pa <= 0.0
        || problem.feed_amounts_mol.size() != component_count
        || !finite(problem.feed_amounts_mol)
        || std::any_of(
            problem.feed_amounts_mol.begin(),
            problem.feed_amounts_mol.end(),
            [](double value) { return value <= 0.0; }
        )) {
        throw std::invalid_argument(
            "phase-equilibrium temperature, pressure, and feed must be positive and finite"
        );
    }
    if ((std::isfinite(transport.source_temperature_min_k)
            && problem.temperature_k < transport.source_temperature_min_k)
        || (std::isfinite(transport.source_temperature_max_k)
            && problem.temperature_k > transport.source_temperature_max_k)) {
        throw std::invalid_argument(
            "phase-equilibrium temperature is outside the model domain"
        );
    }
    if (problem.mode == PhaseEquilibriumMode::Constrained
        && (problem.constrained_roles.empty()
            || problem.constrained_roles.size() > 3)) {
        throw std::invalid_argument(
            "constrained phase equilibrium requires one to three declared roles"
        );
    }
    if (problem.mode == PhaseEquilibriumMode::Discovery
        && !problem.constrained_roles.empty()) {
        throw std::invalid_argument(
            "discovery phase equilibrium cannot contain constrained roles"
        );
    }
    double charge = 0.0;
    double feed_scale = 0.0;
    double ionic_feed = 0.0;
    for (std::size_t component = 0; component < component_count; ++component) {
        charge += charges[component] * problem.feed_amounts_mol[component];
        feed_scale += problem.feed_amounts_mol[component];
        if (charges[component] != 0.0) {
            ionic_feed += problem.feed_amounts_mol[component];
        }
    }
    if (std::abs(charge) > kChargeToleranceMol * std::max(1.0, feed_scale)) {
        throw std::invalid_argument("charged phase-equilibrium feed is not electroneutral");
    }
    if (std::isfinite(transport.total_ion_mole_fraction_max)
        && ionic_feed / feed_scale
            > transport.total_ion_mole_fraction_max + 1.0e-12) {
        throw std::invalid_argument(
            "phase-equilibrium feed ion fraction is outside the model domain"
        );
    }
    for (const auto& hint : problem.hints) {
        if (hint.mole_fractions.size() != component_count
            || !finite(hint.mole_fractions)
            || std::any_of(
                hint.mole_fractions.begin(),
                hint.mole_fractions.end(),
                [](double value) { return value <= 0.0; }
            )
            || !std::isfinite(hint.molar_density_mol_m3)
            || hint.molar_density_mol_m3 <= 0.0
            || std::abs(std::accumulate(
                hint.mole_fractions.begin(), hint.mole_fractions.end(), 0.0
            ) - 1.0) > 1.0e-10) {
            throw std::invalid_argument("phase hint is not a normalized interior state");
        }
        const double packing_fraction = provider.evaluate_packing_fraction(
            problem.temperature_k,
            hint.mole_fractions,
            1.0 / hint.molar_density_mol_m3
        ).value;
        if (!std::isfinite(packing_fraction)
            || packing_fraction <= kPackingFractionMinimum
            || packing_fraction >= kPackingFractionMaximum) {
            throw std::invalid_argument("phase hint density is outside the model domain");
        }
    }
}

}  // namespace

const char* phase_role_name(PhaseRole role) noexcept {
    return role == PhaseRole::Liquid ? "liquid" : "vapor";
}

const char* phase_mode_name(PhaseEquilibriumMode mode) noexcept {
    return mode == PhaseEquilibriumMode::Discovery ? "discovery" : "constrained";
}

PhaseEquilibriumResult solve_phase_equilibrium(
    const ProviderContext& provider,
    const PhaseEquilibriumProblem& problem
) {
    PhaseEquilibriumResult result;
    result.problem = problem;
    result.parameter_fingerprint = provider.fingerprint();
    result.diagnostics = base_diagnostics();
    const auto ids = component_ids(provider);
    const auto charges = component_charges(provider);
    try {
        validate_problem(provider, problem, charges);
    } catch (const std::exception& error) {
        result.diagnostics.solver_status = "input_rejected";
        result.diagnostics.numerical_status = "not_evaluated";
        result.diagnostics.physical_status = "rejected";
        result.diagnostics.failure_stage = "input";
        result.diagnostics.failure_reason = error.what();
        return result;
    }

    const auto topologies = candidate_topologies(problem);
    result.diagnostics.topology_budget = static_cast<int>(topologies.size());
    result.diagnostics.start_budget_per_topology = kStartBudgetPerTopology;
    result.diagnostics.hint_count = static_cast<int>(problem.hints.size());
    result.diagnostics.excluded_topologies = {
        "more_than_one_vapor_phase",
        "solid_phases",
        "reactive_multiphase_equilibrium",
        "more_than_three_phases",
    };
    std::optional<Candidate> best_certified;
    std::optional<Candidate> best_finite;
    std::string first_provider_failure;
    for (const auto& roles : topologies) {
        result.diagnostics.searched_topologies.push_back(topology_name(roles));
        ++result.diagnostics.topology_count;
        auto starts = make_starts(problem.feed_amounts_mol, charges, roles);
        result.diagnostics.hint_reordered_start_count += reorder_starts_from_hints(
            starts,
            provider,
            problem,
            roles,
            problem.hints,
            charges
        );
        for (const auto& start : starts) {
            ++result.diagnostics.attempted_start_count;
            Attempt attempt = run_ipopt(provider, problem, roles, charges, start);
            result.diagnostics.optimizer_iteration_count += attempt.iterations;
            if (first_provider_failure.empty() && !attempt.callback_error.empty()) {
                first_provider_failure = attempt.callback_error;
            }
            if (attempt.variables.empty() || !finite(attempt.variables)) {
                continue;
            }
            try {
                canonicalize_attempt(attempt, problem, charges, provider);
                ProblemEvaluation evaluation = evaluate_problem(
                    provider, problem, attempt.roles, charges, attempt.variables
                );
                PhaseEquilibriumDiagnostics diagnostics = certify_attempt(
                    problem, charges, attempt, evaluation, ids
                );
                Candidate candidate{
                    std::move(attempt),
                    std::move(evaluation),
                    std::move(diagnostics),
                    {},
                    false,
                };
                candidate.states = states_from_evaluation(
                    candidate.evaluation, candidate.attempt.roles
                );
                candidate.certified = candidate.diagnostics.failure_reason.empty();
                if (!best_finite
                    || candidate.evaluation.objective < best_finite->evaluation.objective) {
                    best_finite = candidate;
                }
                if (candidate.certified) {
                    ++result.diagnostics.certified_candidate_count;
                    if (!best_certified
                        || candidate.evaluation.objective
                            < best_certified->evaluation.objective) {
                        best_certified = std::move(candidate);
                    }
                }
            } catch (const std::exception&) {
                continue;
            }
        }
    }
    result.diagnostics.completeness_criterion =
        "all declared role topologies and deterministic starts were evaluated within the recorded finite budgets";
    if (!best_certified) {
        result.diagnostics.solver_status = "search_exhausted";
        result.diagnostics.numerical_status = "uncertified";
        result.diagnostics.physical_status = "rejected";
        result.diagnostics.failure_stage = "certificate";
        result.diagnostics.failure_reason = "no_certified_candidate_in_finite_search";
        result.diagnostics.failure_detail = first_provider_failure;
        if (best_finite) {
            result.best_candidate = best_finite->states;
            const auto retained = result.diagnostics;
            result.diagnostics = best_finite->diagnostics;
            result.diagnostics.topology_budget = retained.topology_budget;
            result.diagnostics.topology_count = retained.topology_count;
            result.diagnostics.start_budget_per_topology =
                retained.start_budget_per_topology;
            result.diagnostics.attempted_start_count = retained.attempted_start_count;
            result.diagnostics.certified_candidate_count = 0;
            result.diagnostics.optimizer_iteration_count =
                retained.optimizer_iteration_count;
            result.diagnostics.hint_count = retained.hint_count;
            result.diagnostics.hint_reordered_start_count =
                retained.hint_reordered_start_count;
            result.diagnostics.searched_topologies = retained.searched_topologies;
            result.diagnostics.excluded_topologies = retained.excluded_topologies;
            result.diagnostics.completeness_criterion = retained.completeness_criterion;
            result.diagnostics.solver_status = "search_exhausted";
            result.diagnostics.numerical_status = "uncertified";
            result.diagnostics.physical_status = "rejected";
            result.diagnostics.failure_stage = "certificate";
            result.diagnostics.failure_reason =
                "no_certified_candidate_in_finite_search";
            result.diagnostics.failure_detail = first_provider_failure;
        }
        return result;
    }

    const auto retained = result.diagnostics;
    result.diagnostics = best_certified->diagnostics;
    result.diagnostics.topology_budget = retained.topology_budget;
    result.diagnostics.topology_count = retained.topology_count;
    result.diagnostics.start_budget_per_topology = retained.start_budget_per_topology;
    result.diagnostics.attempted_start_count = retained.attempted_start_count;
    result.diagnostics.certified_candidate_count = retained.certified_candidate_count;
    result.diagnostics.optimizer_iteration_count = retained.optimizer_iteration_count;
    result.diagnostics.hint_count = retained.hint_count;
    result.diagnostics.hint_reordered_start_count = retained.hint_reordered_start_count;
    result.diagnostics.searched_topologies = retained.searched_topologies;
    result.diagnostics.excluded_topologies = retained.excluded_topologies;
    result.diagnostics.completeness_criterion = retained.completeness_criterion;
    result.diagnostics.certificate_scope = problem.mode
        == PhaseEquilibriumMode::Constrained
        ? "declared_topology_only"
        : "finite_topology_and_start_search";
    result.diagnostics.phase_selection_claim = problem.mode
        == PhaseEquilibriumMode::Constrained
        ? "locally_certified_declared_topology"
        : "lowest_certified_local_candidate_in_recorded_finite_search";
    result.diagnostics.globality_claim = "not_established";
    result.phases = best_certified->states;
    result.total_free_energy_over_rt = best_certified->evaluation.objective;
    if (!provider.active_parameters().identities.empty()) {
        result.sensitivity = evaluate_kkt_sensitivity(
            provider, problem, charges, *best_certified, ids
        );
    }
    return result;
}

}  // namespace epcsaft_engine::equilibrium
