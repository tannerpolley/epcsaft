#include "phase_equilibrium.hpp"

#include "held2.hpp"
#include "held2_tolerances.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Dense>

namespace epcsaft_engine::equilibrium {
namespace {

constexpr double kGasConstant = 8.31446261815324;
constexpr double kPotentialTolerance = 2.0e-8;
constexpr double kPressureTolerance = 1.0e-8;
constexpr double kDistinctTolerance = 1.0e-5;
constexpr int kIterationLimit = 100;

std::vector<double> chart_coordinates(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical
) {
    const auto modified = held2_transform_physical_fractions(
        coordinates, physical
    );
    std::vector<double> result;
    for (std::size_t eos : coordinates.independent_indices) {
        const auto found = std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            eos
        );
        result.push_back(modified[static_cast<std::size_t>(
            found - coordinates.retained_indices.begin()
        )]);
    }
    return result;
}

struct LocalPhase {
    std::vector<double> composition;
    double volume = 0.0;
    epcsaft::NativePhaseState block;
};

LocalPhase evaluate_local(
    const epcsaft::NativeModel& eos,
    const Held2Coordinates& coordinates,
    double temperature,
    const std::vector<double>& independent,
    double log_volume
) {
    LocalPhase result;
    result.composition = held2_lift_independent_fractions(
        coordinates, independent
    );
    result.volume = std::exp(log_volume);
    result.block = eos.evaluate_phase_state(
        temperature, result.composition, result.volume
    );
    return result;
}

PhaseEquilibriumState make_state(
    const epcsaft::NativeModel& eos,
    PhaseRole role,
    PhaseSupport support,
    const LocalPhase& local
) {
    return make_phase_equilibrium_state(
        eos,
        role,
        support,
        local.composition,
        local.volume,
        local.block
    );
}

double volume_seed(
    PhaseRole role,
    const std::array<double, 2>& bounds
) {
    const double lower = std::log(bounds[0]);
    const double upper = std::log(bounds[1]);
    if (role == PhaseRole::Liquid) {
        return std::min(bounds[0] * 1.35, std::nextafter(bounds[1], bounds[0]));
    }
    return std::exp(lower + 0.92 * (upper - lower));
}

std::vector<double> projected_composition(
    const std::vector<double>& reference,
    const std::vector<std::size_t>& mapping
) {
    std::vector<double> result;
    result.reserve(mapping.size());
    for (std::size_t index : mapping) {
        result.push_back(reference[index]);
    }
    double total = std::accumulate(result.begin(), result.end(), 0.0);
    if (!(total > 0.0)) {
        result.assign(result.size(), 1.0 / static_cast<double>(result.size()));
        return result;
    }
    for (double& value : result) {
        value = std::max(value / total, kHeld2ModifiedLowerScale);
    }
    total = std::accumulate(result.begin(), result.end(), 0.0);
    for (double& value : result) {
        value /= total;
    }
    return result;
}

struct Evaluation {
    LocalPhase reference;
    LocalPhase partner;
    std::vector<double> residual;
    Eigen::MatrixXd jacobian;
    std::size_t eos_evaluation_count = 0;
};

struct PotentialRows {
    std::size_t count = 0;
    std::vector<double> reference;
    std::vector<double> partner;
};

std::vector<double> modified_potential_transform(
    const Held2Coordinates& coordinates
) {
    const std::size_t component_count = coordinates.charges.size();
    if (coordinates.eliminated_index == component_count) {
        std::vector<double> result(
            component_count * component_count, 0.0
        );
        for (std::size_t component = 0;
             component < component_count; ++component) {
            result[component * component_count + component] = 1.0;
        }
        return result;
    }
    const std::size_t row_count = coordinates.retained_indices.size();
    std::vector<double> result(row_count * component_count, 0.0);
    const std::size_t eliminated = coordinates.eliminated_index;
    const double eliminated_charge = coordinates.charges[eliminated];
    for (std::size_t row = 0; row < row_count; ++row) {
        const std::size_t component = coordinates.retained_indices[row];
        const double factor = coordinates.modified_factors[row];
        result[row * component_count + component] = 1.0 / factor;
        result[row * component_count + eliminated] =
            -coordinates.charges[component]
            / (eliminated_charge * factor);
    }
    return result;
}

PotentialRows potential_rows(
    const Held2Coordinates& reference_coordinates,
    const Held2Coordinates& partner_coordinates,
    const std::vector<std::size_t>& partner_to_reference
) {
    PotentialRows result;
    const std::size_t reference_count = reference_coordinates.charges.size();
    const std::size_t partner_count = partner_coordinates.charges.size();
    if (phase_topology_is_charged(partner_coordinates.charges)) {
        result.count = partner_coordinates.retained_indices.size();
        result.reference = modified_potential_transform(
            reference_coordinates
        );
        result.partner = modified_potential_transform(partner_coordinates);
        if (result.reference.size() != result.count * reference_count
            || result.partner.size() != result.count * partner_count) {
            throw std::invalid_argument(
                "charged coexistence modified-potential layouts differ"
            );
        }
        return result;
    }
    result.count = partner_count;
    result.reference.assign(result.count * reference_count, 0.0);
    result.partner.assign(result.count * partner_count, 0.0);
    for (std::size_t row = 0; row < result.count; ++row) {
        result.reference[row * reference_count + partner_to_reference[row]] =
            1.0;
        result.partner[row * partner_count + row] = 1.0;
    }
    return result;
}

Evaluation evaluate_system(
    const epcsaft::NativeModel& reference_eos,
    const epcsaft::NativeModel& partner_eos,
    const FixedTCoexistenceProblem& problem,
    const std::vector<double>& reference_charges,
    const std::vector<double>& partner_charges,
    const Held2Coordinates& reference_coordinates,
    const Held2Coordinates& partner_coordinates,
    const std::vector<std::size_t>& partner_to_reference,
    const std::vector<double>& variables
) {
    const std::size_t partner_dimension =
        partner_coordinates.independent_indices.size();
    if (variables.size() != partner_dimension + 2) {
        throw std::invalid_argument("coexistence state dimension changed");
    }
    const std::vector<double> partner_independent(
        variables.begin() + 1,
        variables.begin() + static_cast<std::ptrdiff_t>(
            1 + partner_dimension
        )
    );
    Evaluation result;
    result.reference = evaluate_local(
        reference_eos,
        reference_coordinates,
        problem.temperature_k,
        chart_coordinates(
            reference_coordinates, problem.reference_mole_fractions
        ),
        variables.front()
    );
    result.partner = evaluate_local(
        partner_eos,
        partner_coordinates,
        problem.temperature_k,
        partner_independent,
        variables.back()
    );
    const std::size_t reference_count = reference_charges.size();
    const std::size_t partner_count = partner_charges.size();
    const std::size_t reference_state_count = reference_count + 1;
    const std::size_t partner_state_count = partner_count + 1;
    const PotentialRows potentials = potential_rows(
        reference_coordinates, partner_coordinates, partner_to_reference
    );
    const std::size_t residual_count = potentials.count + 1;
    const auto chart = held2_physical_fraction_jacobian(
        partner_coordinates
    );
    result.residual.assign(residual_count, 0.0);
    result.jacobian = Eigen::MatrixXd::Zero(
        static_cast<Eigen::Index>(residual_count),
        static_cast<Eigen::Index>(variables.size())
    );

    for (std::size_t row = 0; row < potentials.count; ++row) {
        for (std::size_t component = 0;
             component < reference_count; ++component) {
            const double coefficient = potentials.reference[
                row * reference_count + component
            ];
            result.residual[row] += coefficient
                * result.reference.block.gradient[component];
            result.jacobian(
                static_cast<Eigen::Index>(row), 0
            ) += coefficient * result.reference.block.hessian[
                component * reference_state_count + reference_count
            ] * result.reference.volume;
        }
        for (std::size_t component = 0;
             component < partner_count; ++component) {
            result.residual[row] -= potentials.partner[
                row * partner_count + component
            ] * result.partner.block.gradient[component];
        }
        for (std::size_t coordinate = 0;
             coordinate < partner_dimension; ++coordinate) {
            double derivative = 0.0;
            for (std::size_t potential_component = 0;
                 potential_component < partner_count;
                 ++potential_component) {
                for (std::size_t composition_component = 0;
                     composition_component < partner_count;
                     ++composition_component) {
                    derivative += potentials.partner[
                        row * partner_count + potential_component
                    ] * result.partner.block.hessian[
                        potential_component * partner_state_count
                        + composition_component
                    ] * chart[
                        composition_component * partner_dimension + coordinate
                    ];
                }
            }
            result.jacobian(
                static_cast<Eigen::Index>(row),
                static_cast<Eigen::Index>(1 + coordinate)
            ) = -derivative;
        }
        for (std::size_t component = 0;
             component < partner_count; ++component) {
            result.jacobian(
                static_cast<Eigen::Index>(row),
                static_cast<Eigen::Index>(variables.size() - 1)
            ) -= potentials.partner[
                row * partner_count + component
            ] * result.partner.block.hessian[
                component * partner_state_count + partner_count
            ] * result.partner.volume;
        }
    }

    const std::size_t pressure_row = potentials.count;
    const double pressure_scale = std::max({
        1.0,
        std::abs(result.reference.block.pressure_pa),
        std::abs(result.partner.block.pressure_pa),
    });
    result.residual[pressure_row] =
        (result.reference.block.pressure_pa
         - result.partner.block.pressure_pa) / pressure_scale;
    const double rt = kGasConstant * problem.temperature_k;
    result.jacobian(
        static_cast<Eigen::Index>(pressure_row), 0
    ) = -rt * result.reference.block.hessian[
        reference_count * reference_state_count + reference_count
    ] * result.reference.volume / pressure_scale;
    for (std::size_t coordinate = 0;
         coordinate < partner_dimension; ++coordinate) {
        double derivative = 0.0;
        for (std::size_t component = 0;
             component < partner_count; ++component) {
            derivative += result.partner.block.hessian[
                partner_count * partner_state_count + component
            ] * chart[
                component * partner_dimension + coordinate
            ];
        }
        result.jacobian(
            static_cast<Eigen::Index>(pressure_row),
            static_cast<Eigen::Index>(1 + coordinate)
        ) = rt * derivative / pressure_scale;
    }
    result.jacobian(
        static_cast<Eigen::Index>(pressure_row),
        static_cast<Eigen::Index>(variables.size() - 1)
    ) = rt * result.partner.block.hessian[
        partner_count * partner_state_count + partner_count
    ] * result.partner.volume / pressure_scale;
    result.eos_evaluation_count = result.reference.block.evaluations.total()
        + result.partner.block.evaluations.total();
    return result;
}

double norm_inf(const std::vector<double>& values) {
    double result = 0.0;
    for (double value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

bool admissible(
    const epcsaft::NativeModel& reference_eos,
    const epcsaft::NativeModel& partner_eos,
    const FixedTCoexistenceProblem& problem,
    const Held2Coordinates& partner_coordinates,
    const std::vector<double>& variables,
    std::size_t& eos_evaluation_count
) {
    try {
        const std::size_t dimension =
            partner_coordinates.independent_indices.size();
        const std::vector<double> independent(
            variables.begin() + 1,
            variables.begin() + static_cast<std::ptrdiff_t>(1 + dimension)
        );
        const auto composition = held2_lift_independent_fractions(
            partner_coordinates, independent
        );
        const auto reference_bounds =
            reference_eos.evaluate_molar_volume_bounds(
                problem.temperature_k,
                problem.reference_mole_fractions,
                kHeld2PackingFractionMinimum,
                kHeld2PackingFractionMaximum
            );
        const auto partner_bounds =
            partner_eos.evaluate_molar_volume_bounds(
                problem.temperature_k,
                composition,
                kHeld2PackingFractionMinimum,
                kHeld2PackingFractionMaximum
            );
        eos_evaluation_count += reference_bounds.eos_evaluation_count
            + partner_bounds.eos_evaluation_count;
        const double reference_volume = std::exp(variables.front());
        const double partner_volume = std::exp(variables.back());
        return reference_volume > reference_bounds.bounds_m3_per_mol[0]
            && reference_volume < reference_bounds.bounds_m3_per_mol[1]
            && partner_volume > partner_bounds.bounds_m3_per_mol[0]
            && partner_volume < partner_bounds.bounds_m3_per_mol[1];
    } catch (const std::exception&) {
        return false;
    }
}

SolvedStateSensitivity coexistence_sensitivity(
    const epcsaft::NativeModel& reference_eos,
    const epcsaft::NativeModel& partner_eos,
    const FixedTCoexistenceProblem& problem,
    const std::vector<std::size_t>& partner_to_reference,
    const Held2Coordinates& reference_coordinates,
    const Held2Coordinates& partner_coordinates,
    const Evaluation& solved
) {
    SolvedStateSensitivity result;
    const auto reference_active = reference_eos.active_coordinates();
    result.parameter_identities = epcsaft::active_parameter_identities(reference_active);
    result.reference_phase_index = problem.reference_phase_index;
    if (reference_active.empty()) {
        result.failure_reason = "no_active_parameters";
        return result;
    }
    const auto partner_active = partner_eos.active_coordinates();
    const auto reference_identities =
        epcsaft::active_parameter_identities(reference_active);
    const auto partner_identities =
        epcsaft::active_parameter_identities(partner_active);
    const std::size_t parameter_count = reference_active.size();
    const std::size_t partner_count = partner_to_reference.size();
    const std::size_t independent_count =
        partner_coordinates.independent_indices.size();
    result.dimension = independent_count + 2;
    const PotentialRows potentials = potential_rows(
        reference_coordinates, partner_coordinates, partner_to_reference
    );
    if (potentials.count + 1 != result.dimension) {
        result.failure_reason = "coexistence_residual_layout_changed";
        return result;
    }
    const std::size_t reference_count =
        reference_coordinates.charges.size();
    std::vector<std::size_t> partner_parameter_index(
        parameter_count, partner_active.size()
    );
    for (std::size_t parameter = 0;
         parameter < parameter_count; ++parameter) {
        const auto found = std::find(
            partner_identities.begin(), partner_identities.end(),
            reference_identities[parameter]
        );
        if (found != partner_identities.end()) {
            partner_parameter_index[parameter] = static_cast<std::size_t>(
                found - partner_identities.begin()
            );
        }
    }
    epcsaft::NativeParameterizedPhaseState reference;
    std::optional<epcsaft::NativeParameterizedPhaseState> partner;
    try {
        reference = reference_eos.evaluate_parameterized_phase_state(
            problem.temperature_k,
            solved.reference.composition,
            solved.reference.volume,
            reference_active
        );
        result.eos_evaluation_count += reference.evaluations.total();
        if (!partner_active.empty()) {
            partner = partner_eos.evaluate_parameterized_phase_state(
                problem.temperature_k,
                solved.partner.composition,
                solved.partner.volume,
                partner_active
            );
            result.eos_evaluation_count += partner->evaluations.total();
        }
    } catch (const std::exception&) {
        result.failure_reason = "eos_partial_evaluation_failed";
        return result;
    }
    result.residual_state_jacobian.assign(
        result.dimension * result.dimension, 0.0
    );
    for (std::size_t row = 0; row < result.dimension; ++row) {
        for (std::size_t column = 0; column < result.dimension; ++column) {
            result.residual_state_jacobian[
                row * result.dimension + column
            ] = solved.jacobian(
                static_cast<Eigen::Index>(row),
                static_cast<Eigen::Index>(column)
            );
        }
    }
    result.residual_parameter_jacobian.assign(
        result.dimension * parameter_count, 0.0
    );
    const auto partner_ids = phase_eos_component_ids(partner_eos);
    for (std::size_t row = 0; row < potentials.count; ++row) {
        const std::size_t identity_component =
            phase_topology_is_charged(partner_coordinates.charges)
            ? partner_coordinates.retained_indices[row]
            : row;
        result.residual_identities.push_back(
            "phase[" + std::to_string(problem.reference_phase_index)
            + (phase_topology_is_charged(partner_coordinates.charges)
                ? "].modified_potential["
                : "].mu[")
            + partner_ids[identity_component]
            + "]-phase["
            + std::to_string(1 - problem.reference_phase_index)
            + (phase_topology_is_charged(partner_coordinates.charges)
                ? "].modified_potential"
                : "].mu")
        );
        for (std::size_t parameter = 0;
             parameter < parameter_count; ++parameter) {
            double derivative = 0.0;
            for (std::size_t component = 0;
                 component < reference_count; ++component) {
                derivative += potentials.reference[
                    row * reference_count + component
                ] * reference
                    .chemical_potential_parameter_derivatives_over_rt[
                        component * parameter_count + parameter
                    ];
            }
            const std::size_t partner_parameter =
                partner_parameter_index[parameter];
            if (partner.has_value()
                && partner_parameter < partner_active.size()) {
                for (std::size_t component = 0;
                     component < partner_count; ++component) {
                    derivative -= potentials.partner[
                        row * partner_count + component
                    ] * partner
                        ->chemical_potential_parameter_derivatives_over_rt[
                            component * partner_active.size()
                            + partner_parameter
                        ];
                }
            }
            result.residual_parameter_jacobian[
                row * parameter_count + parameter
            ] = derivative;
        }
    }
    result.residual_identities.push_back(
        "phase[" + std::to_string(problem.reference_phase_index)
        + "].pressure_minus_phase["
        + std::to_string(1 - problem.reference_phase_index)
        + "].pressure_scaled"
    );
    const double pressure_scale = std::max({
        1.0,
        std::abs(solved.reference.block.pressure_pa),
        std::abs(solved.partner.block.pressure_pa),
    });
    const std::size_t pressure_row = potentials.count;
    for (std::size_t parameter = 0;
         parameter < parameter_count; ++parameter) {
        double derivative =
            reference.pressure_parameter_derivatives_pa[parameter];
        const std::size_t partner_parameter =
            partner_parameter_index[parameter];
        if (partner.has_value()
            && partner_parameter < partner_active.size()) {
            derivative -= partner->pressure_parameter_derivatives_pa[
                partner_parameter
            ];
        }
        result.residual_parameter_jacobian[
            pressure_row * parameter_count + parameter
        ] = derivative / pressure_scale;
    }
    result.state_identities.push_back(
        "phase[" + std::to_string(problem.reference_phase_index)
        + "].log_molar_volume"
    );
    for (std::size_t coordinate = 0;
         coordinate < independent_count; ++coordinate) {
        result.state_identities.push_back(
            "phase[" + std::to_string(1 - problem.reference_phase_index)
            + "].modified_fraction["
            + partner_ids[partner_coordinates.independent_indices[coordinate]]
            + "]"
        );
    }
    result.state_identities.push_back(
        "phase[" + std::to_string(1 - problem.reference_phase_index)
        + "].log_molar_volume"
    );
    if (!solve_sensitivity_linear_system(result, parameter_count)) {
        return result;
    }
    const std::size_t output_count = 1
        + phase_eos_component_ids(reference_eos).size() + 1
        + partner_ids.size() + 1;
    result.output_parameter_jacobian.assign(
        output_count * parameter_count, 0.0
    );
    std::size_t output = 0;
    result.output_identities.push_back("boundary.pressure_pa");
    const std::size_t reference_state_count = reference_count + 1;
    const double rt = kGasConstant * problem.temperature_k;
    const double pressure_log_volume = -rt
        * reference.phase.hessian[
            reference_count * reference_state_count + reference_count
        ] * solved.reference.volume;
    for (std::size_t parameter = 0;
         parameter < parameter_count; ++parameter) {
        result.output_parameter_jacobian[parameter] =
            reference.pressure_parameter_derivatives_pa[parameter]
            + pressure_log_volume * result.state_parameter_jacobian[
                parameter
            ];
    }
    ++output;
    const auto reference_ids = phase_eos_component_ids(reference_eos);
    for (const auto& id : reference_ids) {
        result.output_identities.push_back(
            "phase[" + std::to_string(problem.reference_phase_index)
            + "].mole_fraction[" + id + "]"
        );
        ++output;
    }
    result.output_identities.push_back(
        "phase[" + std::to_string(problem.reference_phase_index)
        + "].molar_volume_m3_per_mol"
    );
    for (std::size_t parameter = 0;
         parameter < parameter_count; ++parameter) {
        result.output_parameter_jacobian[
            output * parameter_count + parameter
        ] = solved.reference.volume * result.state_parameter_jacobian[
            parameter
        ];
    }
    ++output;
    const auto chart = held2_physical_fraction_jacobian(
        partner_coordinates
    );
    for (std::size_t component = 0;
         component < partner_ids.size(); ++component) {
        result.output_identities.push_back(
            "phase[" + std::to_string(1 - problem.reference_phase_index)
            + "].mole_fraction[" + partner_ids[component] + "]"
        );
        for (std::size_t parameter = 0;
             parameter < parameter_count; ++parameter) {
            double derivative = 0.0;
            for (std::size_t coordinate = 0;
                 coordinate < independent_count; ++coordinate) {
                derivative += chart[
                    component * independent_count + coordinate
                ] * result.state_parameter_jacobian[
                    (1 + coordinate) * parameter_count + parameter
                ];
            }
            result.output_parameter_jacobian[
                output * parameter_count + parameter
            ] = derivative;
        }
        ++output;
    }
    result.output_identities.push_back(
        "phase[" + std::to_string(1 - problem.reference_phase_index)
        + "].molar_volume_m3_per_mol"
    );
    for (std::size_t parameter = 0;
         parameter < parameter_count; ++parameter) {
        result.output_parameter_jacobian[
            output * parameter_count + parameter
        ] = solved.partner.volume * result.state_parameter_jacobian[
            (result.dimension - 1) * parameter_count + parameter
        ];
    }
    result.status = "available";
    result.failure_reason.clear();
    return result;
}

SolvedStateSensitivity coexistence_reference_direction_sensitivity(
    const epcsaft::NativeModel& reference_eos,
    const epcsaft::NativeModel& partner_eos,
    const FixedTCoexistenceProblem& problem,
    const std::vector<std::size_t>& partner_to_reference,
    const Held2Coordinates& reference_coordinates,
    const Held2Coordinates& partner_coordinates,
    const Evaluation& solved
) {
    SolvedStateSensitivity result;
    result.parameter_identities = problem.reference_direction_identities;
    result.reference_phase_index = problem.reference_phase_index;
    const std::size_t direction_count =
        problem.reference_composition_directions.size();
    if (direction_count == 0
        || result.parameter_identities.size() != direction_count) {
        result.failure_reason = "reference_direction_layout_changed";
        return result;
    }
    const std::size_t reference_count =
        reference_coordinates.charges.size();
    const std::size_t independent_count =
        partner_coordinates.independent_indices.size();
    result.dimension = independent_count + 2;
    const PotentialRows potentials = potential_rows(
        reference_coordinates, partner_coordinates, partner_to_reference
    );
    if (potentials.count + 1 != result.dimension) {
        result.failure_reason = "coexistence_residual_layout_changed";
        return result;
    }
    result.residual_state_jacobian.assign(
        result.dimension * result.dimension, 0.0
    );
    for (std::size_t row = 0; row < result.dimension; ++row) {
        for (std::size_t column = 0; column < result.dimension; ++column) {
            result.residual_state_jacobian[
                row * result.dimension + column
            ] = solved.jacobian(
                static_cast<Eigen::Index>(row),
                static_cast<Eigen::Index>(column)
            );
        }
    }
    result.residual_parameter_jacobian.assign(
        result.dimension * direction_count, 0.0
    );
    const auto partner_ids = phase_eos_component_ids(partner_eos);
    const std::size_t reference_state_count = reference_count + 1;
    for (std::size_t row = 0; row < potentials.count; ++row) {
        const std::size_t identity_component =
            phase_topology_is_charged(partner_coordinates.charges)
            ? partner_coordinates.retained_indices[row]
            : row;
        result.residual_identities.push_back(
            "phase[" + std::to_string(problem.reference_phase_index)
            + (phase_topology_is_charged(partner_coordinates.charges)
                ? "].modified_potential["
                : "].mu[")
            + partner_ids[identity_component]
            + "]-phase["
            + std::to_string(1 - problem.reference_phase_index)
            + (phase_topology_is_charged(partner_coordinates.charges)
                ? "].modified_potential"
                : "].mu")
        );
        for (std::size_t direction = 0;
             direction < direction_count; ++direction) {
            double derivative = 0.0;
            for (std::size_t potential_component = 0;
                 potential_component < reference_count;
                 ++potential_component) {
                const double coefficient = potentials.reference[
                    row * reference_count + potential_component
                ];
                for (std::size_t composition_component = 0;
                     composition_component < reference_count;
                     ++composition_component) {
                    derivative += coefficient
                        * solved.reference.block.hessian[
                            potential_component * reference_state_count
                            + composition_component
                        ]
                        * problem.reference_composition_directions[direction][
                            composition_component
                        ];
                }
            }
            result.residual_parameter_jacobian[
                row * direction_count + direction
            ] = derivative;
        }
    }
    result.residual_identities.push_back(
        "phase[" + std::to_string(problem.reference_phase_index)
        + "].pressure_minus_phase["
        + std::to_string(1 - problem.reference_phase_index)
        + "].pressure_scaled"
    );
    const double pressure_scale = std::max({
        1.0,
        std::abs(solved.reference.block.pressure_pa),
        std::abs(solved.partner.block.pressure_pa),
    });
    const double rt = kGasConstant * problem.temperature_k;
    const std::size_t pressure_row = potentials.count;
    for (std::size_t direction = 0;
         direction < direction_count; ++direction) {
        double derivative = 0.0;
        for (std::size_t component = 0;
             component < reference_count; ++component) {
            derivative -= rt * solved.reference.block.hessian[
                reference_count * reference_state_count + component
            ] * problem.reference_composition_directions[direction][component];
        }
        result.residual_parameter_jacobian[
            pressure_row * direction_count + direction
        ] = derivative / pressure_scale;
    }
    result.state_identities.push_back(
        "phase[" + std::to_string(problem.reference_phase_index)
        + "].log_molar_volume"
    );
    for (std::size_t coordinate = 0;
         coordinate < independent_count; ++coordinate) {
        result.state_identities.push_back(
            "phase[" + std::to_string(1 - problem.reference_phase_index)
            + "].modified_fraction["
            + partner_ids[partner_coordinates.independent_indices[coordinate]]
            + "]"
        );
    }
    result.state_identities.push_back(
        "phase[" + std::to_string(1 - problem.reference_phase_index)
        + "].log_molar_volume"
    );
    if (!solve_sensitivity_linear_system(result, direction_count)) {
        return result;
    }
    const auto reference_ids = phase_eos_component_ids(reference_eos);
    const std::size_t output_count = 1 + reference_ids.size() + 1
        + partner_ids.size() + 1;
    result.output_parameter_jacobian.assign(
        output_count * direction_count, 0.0
    );
    result.output_identities.push_back("boundary.pressure_pa");
    const double pressure_log_volume = -rt
        * solved.reference.block.hessian[
            reference_count * reference_state_count + reference_count
        ] * solved.reference.volume;
    for (std::size_t direction = 0;
         direction < direction_count; ++direction) {
        double direct = 0.0;
        for (std::size_t component = 0;
             component < reference_count; ++component) {
            direct -= rt * solved.reference.block.hessian[
                reference_count * reference_state_count + component
            ] * problem.reference_composition_directions[direction][component];
        }
        result.output_parameter_jacobian[direction] = direct
            + pressure_log_volume
                * result.state_parameter_jacobian[direction];
    }
    std::size_t output = 1;
    for (std::size_t component = 0;
         component < reference_ids.size(); ++component) {
        result.output_identities.push_back(
            "phase[" + std::to_string(problem.reference_phase_index)
            + "].mole_fraction[" + reference_ids[component] + "]"
        );
        for (std::size_t direction = 0;
             direction < direction_count; ++direction) {
            result.output_parameter_jacobian[
                output * direction_count + direction
            ] = problem.reference_composition_directions[direction][component];
        }
        ++output;
    }
    result.output_identities.push_back(
        "phase[" + std::to_string(problem.reference_phase_index)
        + "].molar_volume_m3_per_mol"
    );
    for (std::size_t direction = 0;
         direction < direction_count; ++direction) {
        result.output_parameter_jacobian[
            output * direction_count + direction
        ] = solved.reference.volume
            * result.state_parameter_jacobian[direction];
    }
    ++output;
    const auto chart = held2_physical_fraction_jacobian(
        partner_coordinates
    );
    for (std::size_t component = 0;
         component < partner_ids.size(); ++component) {
        result.output_identities.push_back(
            "phase[" + std::to_string(1 - problem.reference_phase_index)
            + "].mole_fraction[" + partner_ids[component] + "]"
        );
        for (std::size_t direction = 0;
             direction < direction_count; ++direction) {
            double derivative = 0.0;
            for (std::size_t coordinate = 0;
                 coordinate < independent_count; ++coordinate) {
                derivative += chart[
                    component * independent_count + coordinate
                ] * result.state_parameter_jacobian[
                    (1 + coordinate) * direction_count + direction
                ];
            }
            result.output_parameter_jacobian[
                output * direction_count + direction
            ] = derivative;
        }
        ++output;
    }
    result.output_identities.push_back(
        "phase[" + std::to_string(1 - problem.reference_phase_index)
        + "].molar_volume_m3_per_mol"
    );
    for (std::size_t direction = 0;
         direction < direction_count; ++direction) {
        result.output_parameter_jacobian[
            output * direction_count + direction
        ] = solved.partner.volume * result.state_parameter_jacobian[
            (result.dimension - 1) * direction_count + direction
        ];
    }
    result.status = "available";
    result.failure_reason.clear();
    return result;
}

SolvedStateSensitivity coexistence_temperature_sensitivity(
    const epcsaft::NativeModel& reference_eos,
    const epcsaft::NativeModel& partner_eos,
    const FixedTCoexistenceProblem& problem,
    const std::vector<std::size_t>& partner_to_reference,
    const Held2Coordinates& reference_coordinates,
    const Held2Coordinates& partner_coordinates,
    const Evaluation& solved
) {
    SolvedStateSensitivity result;
    result.parameter_identities = {"state/temperature_k"};
    result.reference_phase_index = problem.reference_phase_index;
    const std::size_t reference_count = reference_coordinates.charges.size();
    const std::size_t partner_count = partner_to_reference.size();
    const std::size_t independent_count =
        partner_coordinates.independent_indices.size();
    result.dimension = independent_count + 2;
    const PotentialRows potentials = potential_rows(
        reference_coordinates, partner_coordinates, partner_to_reference
    );
    if (potentials.count + 1 != result.dimension) {
        result.failure_reason = "coexistence_residual_layout_changed";
        return result;
    }
    epcsaft::NativeTemperaturePhaseDerivatives reference;
    epcsaft::NativeTemperaturePhaseDerivatives partner;
    try {
        reference = reference_eos.evaluate_temperature_phase(
            problem.temperature_k,
            solved.reference.composition,
            solved.reference.volume
        );
        partner = partner_eos.evaluate_temperature_phase(
            problem.temperature_k,
            solved.partner.composition,
            solved.partner.volume
        );
        result.eos_evaluation_count = reference.evaluations.total()
            + partner.evaluations.total();
    } catch (const std::exception&) {
        result.failure_reason = "eos_temperature_partial_evaluation_failed";
        return result;
    }
    if (reference.chemical_potential_temperature_derivatives_over_rt_per_k.size()
            != reference_count
        || partner.chemical_potential_temperature_derivatives_over_rt_per_k.size()
            != partner_count) {
        result.failure_reason = "temperature_partial_layout_changed";
        return result;
    }
    result.residual_state_jacobian.assign(
        result.dimension * result.dimension, 0.0
    );
    for (std::size_t row = 0; row < result.dimension; ++row) {
        for (std::size_t column = 0; column < result.dimension; ++column) {
            result.residual_state_jacobian[row * result.dimension + column] =
                solved.jacobian(
                    static_cast<Eigen::Index>(row),
                    static_cast<Eigen::Index>(column)
                );
        }
    }
    result.residual_parameter_jacobian.assign(result.dimension, 0.0);
    const auto partner_ids = phase_eos_component_ids(partner_eos);
    for (std::size_t row = 0; row < potentials.count; ++row) {
        const std::size_t identity_component =
            phase_topology_is_charged(partner_coordinates.charges)
            ? partner_coordinates.retained_indices[row]
            : row;
        result.residual_identities.push_back(
            "phase[" + std::to_string(problem.reference_phase_index)
            + (phase_topology_is_charged(partner_coordinates.charges)
                ? "].modified_potential[" : "].mu[")
            + partner_ids[identity_component]
            + "]-phase[" + std::to_string(1 - problem.reference_phase_index)
            + (phase_topology_is_charged(partner_coordinates.charges)
                ? "].modified_potential" : "].mu")
        );
        double derivative = 0.0;
        for (std::size_t component = 0; component < reference_count; ++component) {
            derivative += potentials.reference[row * reference_count + component]
                * reference
                    .chemical_potential_temperature_derivatives_over_rt_per_k[component];
        }
        for (std::size_t component = 0; component < partner_count; ++component) {
            derivative -= potentials.partner[row * partner_count + component]
                * partner
                    .chemical_potential_temperature_derivatives_over_rt_per_k[component];
        }
        result.residual_parameter_jacobian[row] = derivative;
    }
    result.residual_identities.push_back(
        "phase[" + std::to_string(problem.reference_phase_index)
        + "].pressure_minus_phase["
        + std::to_string(1 - problem.reference_phase_index)
        + "].pressure_scaled"
    );
    const double pressure_scale = std::max({
        1.0,
        std::abs(solved.reference.block.pressure_pa),
        std::abs(solved.partner.block.pressure_pa),
    });
    result.residual_parameter_jacobian[potentials.count] =
        (reference.pressure_temperature_derivative_pa_per_k
         - partner.pressure_temperature_derivative_pa_per_k)
        / pressure_scale;
    result.state_identities.push_back(
        "phase[" + std::to_string(problem.reference_phase_index)
        + "].log_molar_volume"
    );
    for (std::size_t coordinate = 0; coordinate < independent_count; ++coordinate) {
        result.state_identities.push_back(
            "phase[" + std::to_string(1 - problem.reference_phase_index)
            + "].modified_fraction["
            + partner_ids[partner_coordinates.independent_indices[coordinate]]
            + "]"
        );
    }
    result.state_identities.push_back(
        "phase[" + std::to_string(1 - problem.reference_phase_index)
        + "].log_molar_volume"
    );
    if (!solve_sensitivity_linear_system(result, 1)) {
        return result;
    }
    const auto reference_ids = phase_eos_component_ids(reference_eos);
    const std::size_t output_count = 1 + reference_ids.size() + 1
        + partner_ids.size() + 1;
    result.output_parameter_jacobian.assign(output_count, 0.0);
    std::size_t output = 0;
    result.output_identities.push_back("boundary.pressure_pa");
    const std::size_t reference_state_count = reference_count + 1;
    const double rt = kGasConstant * problem.temperature_k;
    const double pressure_log_volume = -rt
        * solved.reference.block.hessian[
            reference_count * reference_state_count + reference_count
        ] * solved.reference.volume;
    result.output_parameter_jacobian[output++] =
        reference.pressure_temperature_derivative_pa_per_k
        + pressure_log_volume * result.state_parameter_jacobian[0];
    for (const auto& id : reference_ids) {
        result.output_identities.push_back(
            "phase[" + std::to_string(problem.reference_phase_index)
            + "].mole_fraction[" + id + "]"
        );
        ++output;
    }
    result.output_identities.push_back(
        "phase[" + std::to_string(problem.reference_phase_index)
        + "].molar_volume_m3_per_mol"
    );
    result.output_parameter_jacobian[output++] =
        solved.reference.volume * result.state_parameter_jacobian[0];
    const auto chart = held2_physical_fraction_jacobian(partner_coordinates);
    for (std::size_t component = 0; component < partner_ids.size(); ++component) {
        result.output_identities.push_back(
            "phase[" + std::to_string(1 - problem.reference_phase_index)
            + "].mole_fraction[" + partner_ids[component] + "]"
        );
        for (std::size_t coordinate = 0; coordinate < independent_count; ++coordinate) {
            result.output_parameter_jacobian[output] +=
                chart[component * independent_count + coordinate]
                * result.state_parameter_jacobian[1 + coordinate];
        }
        ++output;
    }
    result.output_identities.push_back(
        "phase[" + std::to_string(1 - problem.reference_phase_index)
        + "].molar_volume_m3_per_mol"
    );
    result.output_parameter_jacobian[output] = solved.partner.volume
        * result.state_parameter_jacobian[result.dimension - 1];
    result.status = "available";
    result.failure_reason.clear();
    return result;
}

void validate_problem(
    const epcsaft::NativeModel& reference_eos,
    const epcsaft::NativeModel& partner_eos,
    const FixedTCoexistenceProblem& problem,
    std::vector<std::size_t>& partner_to_reference
) {
    const auto reference_ids = phase_eos_component_ids(reference_eos);
    const auto partner_ids = phase_eos_component_ids(partner_eos);
    const auto reference_charges = phase_eos_component_charges(
        reference_eos
    );
    const auto partner_charges = phase_eos_component_charges(
        partner_eos
    );
    if (!std::isfinite(problem.temperature_k)
        || problem.temperature_k <= 0.0
        || problem.phases.size() != 2
        || problem.reference_phase_index >= 2
        || problem.reference_mole_fractions.size() != reference_ids.size()
        || std::any_of(
            problem.reference_mole_fractions.begin(),
            problem.reference_mole_fractions.end(),
            [](double value) {
                return !std::isfinite(value) || value <= 0.0;
            }
        )
        || std::abs(
            std::accumulate(
                problem.reference_mole_fractions.begin(),
                problem.reference_mole_fractions.end(),
                0.0
            ) - 1.0
        ) > 1.0e-12
        || problem.phases[problem.reference_phase_index].support
            != PhaseSupport::AllComponents) {
        throw std::invalid_argument(
            "fixed-T coexistence requires two phases and a normalized all-component reference phase"
        );
    }
    if (problem.reference_direction_identities.size()
            != problem.reference_composition_directions.size()
        || std::any_of(
            problem.reference_direction_identities.begin(),
            problem.reference_direction_identities.end(),
            [](const auto& identity) { return identity.empty(); }
        )) {
        throw std::invalid_argument(
            "reference-composition direction identities and rows must match"
        );
    }
    {
        auto identities = problem.reference_direction_identities;
        std::sort(identities.begin(), identities.end());
        if (std::adjacent_find(identities.begin(), identities.end())
            != identities.end()) {
            throw std::invalid_argument(
                "reference-composition direction identities must be unique"
            );
        }
    }
    for (const auto& direction : problem.reference_composition_directions) {
        if (direction.size() != reference_ids.size()
            || std::any_of(
                direction.begin(), direction.end(),
                [](double value) { return !std::isfinite(value); }
            )) {
            throw std::invalid_argument(
                "reference-composition directions must be finite full-component rows"
            );
        }
        const double scale = std::max(
            1.0,
            std::accumulate(
                direction.begin(), direction.end(), 0.0,
                [](double total, double value) {
                    return total + std::abs(value);
                }
            )
        );
        const double amount_tangent = std::accumulate(
            direction.begin(), direction.end(), 0.0
        );
        double charge_tangent = 0.0;
        for (std::size_t component = 0;
             component < direction.size(); ++component) {
            charge_tangent += direction[component]
                * reference_charges[component];
        }
        if (std::abs(amount_tangent) > 1.0e-10 * scale
            || std::abs(charge_tangent) > 1.0e-10 * scale) {
            throw std::invalid_argument(
                "reference-composition directions must preserve normalization and electroneutrality"
            );
        }
    }
    const std::size_t partner_phase = 1 - problem.reference_phase_index;
    const PhaseSupport support = problem.phases[partner_phase].support;
    partner_to_reference.clear();
    for (const auto& id : partner_ids) {
        const auto found = std::find(reference_ids.begin(), reference_ids.end(), id);
        if (found == reference_ids.end()) {
            throw std::invalid_argument(
                "partner component is absent from the reference model"
            );
        }
        partner_to_reference.push_back(static_cast<std::size_t>(
            found - reference_ids.begin()
        ));
    }
    if (support == PhaseSupport::AllComponents && partner_ids != reference_ids) {
        throw std::invalid_argument(
            "all-component coexistence requires identical EOS component order"
        );
    }
    if (support == PhaseSupport::AllComponents
        && partner_charges != reference_charges) {
        throw std::invalid_argument(
            "all-component coexistence requires identical charge topology"
        );
    }
    if (support == PhaseSupport::NeutralComponents) {
        std::vector<std::string> expected;
        for (std::size_t component = 0;
             component < reference_ids.size(); ++component) {
            if (reference_charges[component] == 0.0) {
                expected.push_back(reference_ids[component]);
            }
        }
        if (partner_ids != expected) {
            throw std::invalid_argument(
                "neutral support must contain exactly the reference model neutral components"
            );
        }
    }
    if (support == PhaseSupport::DeclaredNeutralComponents
        && (partner_ids.empty()
            || std::any_of(
                partner_charges.begin(), partner_charges.end(),
                [](double charge) { return charge != 0.0; }
            ))) {
        throw std::invalid_argument(
            "declared-neutral support must contain only neutral reference components"
        );
    }
    double charge = 0.0;
    double scale = 0.0;
    for (std::size_t component = 0;
         component < reference_charges.size(); ++component) {
        charge += reference_charges[component]
            * problem.reference_mole_fractions[component];
        scale += std::abs(reference_charges[component])
            * problem.reference_mole_fractions[component];
    }
    if (!audit_held2_tolerance(kHeld2ChargeBalance, charge, scale).passed) {
        throw std::invalid_argument(
            "coexistence reference composition is not electroneutral"
        );
    }
    std::vector<bool> hinted(2, false);
    for (const auto& hint : problem.hints) {
        if (hint.phase_index >= 2 || hinted[hint.phase_index]
            || hint.mole_fractions.size()
                != (hint.phase_index == problem.reference_phase_index
                    ? reference_ids.size()
                    : partner_ids.size())
            || !std::isfinite(hint.molar_density_mol_m3)
            || hint.molar_density_mol_m3 <= 0.0
            || std::any_of(
                hint.mole_fractions.begin(),
                hint.mole_fractions.end(),
                [](double value) {
                    return !std::isfinite(value) || value <= 0.0;
                }
            )
            || std::abs(
                std::accumulate(
                    hint.mole_fractions.begin(),
                    hint.mole_fractions.end(),
                    0.0
                ) - 1.0
            ) > 1.0e-12) {
            throw std::invalid_argument(
                "coexistence phase hints must be unique, normalized, positive, and finite"
            );
        }
        hinted[hint.phase_index] = true;
    }
}

}  // namespace

FixedTCoexistenceResult solve_fixed_t_coexistence(
    const epcsaft::NativeModel& reference_eos,
    const epcsaft::NativeModel& partner_eos,
    const FixedTCoexistenceProblem& problem
) {
    FixedTCoexistenceResult result;
    result.reference_phase_index = problem.reference_phase_index;
    result.parameter_fingerprint = reference_eos.parameter_fingerprint();
    result.partner_parameter_fingerprint = partner_eos.parameter_fingerprint();
    result.partner_component_ids = phase_eos_component_ids(partner_eos);
    result.diagnostics.route = "held_fixed_t_coexistence";
    result.diagnostics.certificate_scope =
        "fixed_t_declared_coexistence";
    result.diagnostics.phase_selection_claim =
        "declared_support_constrained_coexistence";
    result.diagnostics.search_scope =
        "fixed_reference_composition_neutral_enrichment_matrix";
    result.diagnostics.search_budget =
        "enrichment_levels=0,0.35,0.60,0.85;newton_iterations_per_start=100";
    result.diagnostics.search_exclusions = {
        "undeclared_partner_support",
        "additional_phase_counts",
    };
    result.diagnostics.completeness_criterion =
        "unique_certified_root_among_finite_deterministic_starts;globality_not_established";
    result.diagnostics.role_evidence = {
        "native_pressure_derivative_strict_stability",
        "ordered_molar_volume_branch",
    };
    result.diagnostics.hint_count = static_cast<int>(problem.hints.size());
    std::size_t eos_evaluation_count = 0;
    try {
        std::vector<std::size_t> partner_to_reference;
        validate_problem(
            reference_eos,
            partner_eos,
            problem,
            partner_to_reference
        );
        const auto reference_charges = phase_eos_component_charges(
            reference_eos
        );
        const auto partner_charges = phase_eos_component_charges(
            partner_eos
        );
        const auto reference_coordinates = make_held2_coordinates(
            reference_charges, phase_eos_component_ids(reference_eos)
        );
        const auto partner_coordinates = make_held2_coordinates(
            partner_charges, phase_eos_component_ids(partner_eos)
        );
        const auto partner_composition = projected_composition(
            problem.reference_mole_fractions, partner_to_reference
        );
        const auto reference_bounds =
            reference_eos.evaluate_molar_volume_bounds(
                problem.temperature_k,
                problem.reference_mole_fractions,
                kHeld2PackingFractionMinimum,
                kHeld2PackingFractionMaximum
            );
        eos_evaluation_count += reference_bounds.eos_evaluation_count;
        const std::size_t partner_phase = 1 - problem.reference_phase_index;
        const double reference_volume = volume_seed(
            problem.phases[problem.reference_phase_index].role,
            reference_bounds.bounds_m3_per_mol
        );
        std::vector<std::vector<double>> start_compositions{
            partner_composition,
        };
        for (std::size_t component = 0;
             component < partner_composition.size(); ++component) {
            if (partner_charges[component] != 0.0) {
                continue;
            }
            for (double enrichment : {0.35, 0.60, 0.85}) {
                std::vector<double> enriched = partner_composition;
                for (double& value : enriched) {
                    value *= 1.0 - enrichment;
                }
                enriched[component] += enrichment;
                const bool duplicate = std::any_of(
                    start_compositions.begin(),
                    start_compositions.end(),
                    [&enriched](const auto& existing) {
                        double difference = 0.0;
                        for (std::size_t index = 0;
                             index < enriched.size(); ++index) {
                            difference = std::max(
                                difference,
                                std::abs(enriched[index] - existing[index])
                            );
                        }
                        return difference <= 1.0e-12;
                    }
                );
                if (!duplicate) {
                    start_compositions.push_back(std::move(enriched));
                }
            }
        }
        const PhaseHint* partner_hint = nullptr;
        for (const auto& hint : problem.hints) {
            if (hint.phase_index == partner_phase) {
                partner_hint = &hint;
            }
        }
        if (partner_hint != nullptr) {
            const auto target = partner_hint->mole_fractions;
            const auto before = start_compositions;
            std::stable_sort(
                start_compositions.begin(),
                start_compositions.end(),
                [&target](const auto& left, const auto& right) {
                    double left_distance = 0.0;
                    double right_distance = 0.0;
                    for (std::size_t component = 0;
                         component < target.size(); ++component) {
                        left_distance += std::abs(
                            left[component] - target[component]
                        );
                        right_distance += std::abs(
                            right[component] - target[component]
                        );
                    }
                    return left_distance < right_distance;
                }
            );
            result.diagnostics.hint_reordered_start_count =
                start_compositions == before ? 0 : 1;
        }
        const auto solve_start = [&](std::vector<double> variables)
            -> std::optional<Evaluation> {
            ++result.diagnostics.optimizer_solve_count;
            try {
                for (int iteration = 0;
                     iteration < kIterationLimit; ++iteration) {
                    ++result.diagnostics.optimizer_iteration_count;
                    Evaluation evaluation = evaluate_system(
                        reference_eos,
                        partner_eos,
                        problem,
                        reference_charges,
                        partner_charges,
                        reference_coordinates,
                        partner_coordinates,
                        partner_to_reference,
                        variables
                    );
                    eos_evaluation_count += evaluation.eos_evaluation_count;
                    const double current_norm = norm_inf(
                        evaluation.residual
                    );
                    double potential_norm = 0.0;
                    for (std::size_t row = 0;
                         row + 1 < evaluation.residual.size(); ++row) {
                        potential_norm = std::max(
                            potential_norm,
                            std::abs(evaluation.residual[row])
                        );
                    }
                    if (potential_norm <= kPotentialTolerance
                        && std::abs(evaluation.residual.back())
                            <= kPressureTolerance) {
                        return evaluation;
                    }
                    Eigen::VectorXd residual(
                        static_cast<Eigen::Index>(
                            evaluation.residual.size()
                        )
                    );
                    for (std::size_t row = 0;
                         row < evaluation.residual.size(); ++row) {
                        residual[static_cast<Eigen::Index>(row)] =
                            evaluation.residual[row];
                    }
                    Eigen::FullPivLU<Eigen::MatrixXd> decomposition(
                        evaluation.jacobian
                    );
                    if (decomposition.rank()
                        != evaluation.jacobian.rows()) {
                        return std::nullopt;
                    }
                    const Eigen::VectorXd step =
                        decomposition.solve(-residual);
                    bool improved = false;
                    double length = 1.0;
                    for (int backtrack = 0;
                         backtrack < 30; ++backtrack) {
                        std::vector<double> trial = variables;
                        for (std::size_t coordinate = 0;
                             coordinate < trial.size(); ++coordinate) {
                            trial[coordinate] += length * step[
                                static_cast<Eigen::Index>(coordinate)
                            ];
                        }
                        if (admissible(
                                reference_eos,
                                partner_eos,
                                problem,
                                partner_coordinates,
                                trial,
                                eos_evaluation_count
                            )) {
                            try {
                                const auto candidate = evaluate_system(
                                    reference_eos,
                                    partner_eos,
                                    problem,
                                    reference_charges,
                                    partner_charges,
                                    reference_coordinates,
                                    partner_coordinates,
                                    partner_to_reference,
                                    trial
                                );
                                eos_evaluation_count += candidate.eos_evaluation_count;
                                if (norm_inf(candidate.residual)
                                    < current_norm) {
                                    variables = std::move(trial);
                                    improved = true;
                                    break;
                                }
                            } catch (const std::exception&) {
                            }
                        }
                        length *= 0.5;
                    }
                    if (!improved) {
                        return std::nullopt;
                    }
                }
            } catch (const std::exception&) {
            }
            return std::nullopt;
        };
        std::vector<Evaluation> accepted;
        const auto same_root = [](const Evaluation& left,
                                  const Evaluation& right) {
            const double left_pressure = 0.5 * (
                left.reference.block.pressure_pa
                + left.partner.block.pressure_pa
            );
            const double right_pressure = 0.5 * (
                right.reference.block.pressure_pa
                + right.partner.block.pressure_pa
            );
            if (std::abs(std::log(left_pressure / right_pressure)) > 2.0e-8
                || std::abs(std::log(
                    left.reference.volume / right.reference.volume
                )) > 2.0e-8
                || std::abs(std::log(
                    left.partner.volume / right.partner.volume
                )) > 2.0e-8) {
                return false;
            }
            for (std::size_t component = 0;
                 component < left.partner.composition.size(); ++component) {
                if (std::abs(
                        left.partner.composition[component]
                        - right.partner.composition[component]
                    ) > 2.0e-8) {
                    return false;
                }
            }
            return true;
        };
        for (const auto& composition : start_compositions) {
            const auto chart = chart_coordinates(
                partner_coordinates, composition
            );
            const auto bounds = partner_eos.evaluate_molar_volume_bounds(
                problem.temperature_k,
                composition,
                kHeld2PackingFractionMinimum,
                kHeld2PackingFractionMaximum
            );
            eos_evaluation_count += bounds.eos_evaluation_count;
            std::vector<double> variables{
                std::log(reference_volume),
            };
            variables.insert(variables.end(), chart.begin(), chart.end());
            const double partner_volume =
                problem.phases[partner_phase].role == PhaseRole::Liquid
                ? std::exp(
                    std::log(bounds.bounds_m3_per_mol[0])
                    + 0.05 * (
                        std::log(bounds.bounds_m3_per_mol[1])
                        - std::log(bounds.bounds_m3_per_mol[0])
                    )
                )
                : volume_seed(
                    problem.phases[partner_phase].role,
                    bounds.bounds_m3_per_mol
                );
            variables.push_back(std::log(partner_volume));
            auto candidate = solve_start(std::move(variables));
            if (!candidate.has_value()) {
                continue;
            }
            const auto reference_state = make_state(
                reference_eos,
                problem.phases[problem.reference_phase_index].role,
                PhaseSupport::AllComponents,
                candidate->reference
            );
            const auto partner_state = make_state(
                partner_eos,
                problem.phases[partner_phase].role,
                problem.phases[partner_phase].support,
                candidate->partner
            );
            const bool stable =
                reference_state.mechanical_class == "strict_stable"
                && partner_state.mechanical_class == "strict_stable";
            const bool ordered = reference_state.role == partner_state.role
                || (reference_state.role == PhaseRole::Liquid
                    ? reference_state.molar_volume_m3_per_mol
                        < partner_state.molar_volume_m3_per_mol
                    : reference_state.molar_volume_m3_per_mol
                        > partner_state.molar_volume_m3_per_mol);
            const double distinctness = std::abs(std::log(
                reference_state.molar_volume_m3_per_mol
                / partner_state.molar_volume_m3_per_mol
            ));
            if (stable && ordered && distinctness > kDistinctTolerance
                && candidate->reference.block.pressure_pa > 0.0
                && candidate->partner.block.pressure_pa > 0.0) {
                ++result.diagnostics.accepted_start_count;
                if (!std::any_of(
                        accepted.begin(), accepted.end(),
                        [&candidate, &same_root](const auto& existing) {
                            return same_root(existing, *candidate);
                        }
                    )) {
                    accepted.push_back(std::move(*candidate));
                }
            }
        }
        result.diagnostics.distinct_root_count = static_cast<int>(
            accepted.size()
        );
        if (accepted.empty()) {
            throw std::runtime_error(
                "no deterministic coexistence start produced a certified branch"
            );
        }
        if (accepted.size() != 1) {
            result.diagnostics.eos_evaluation_count =
                static_cast<int>(eos_evaluation_count);
            result.diagnostics.solver_status = "converged";
            result.diagnostics.numerical_status = "certified";
            result.diagnostics.physical_status = "rejected";
            result.diagnostics.failure_stage = "branch_selection";
            result.diagnostics.failure_reason =
                "multiple_unresolved_coexistence_branches";
            return result;
        }
        result.diagnostics.eos_evaluation_count =
            static_cast<int>(eos_evaluation_count);
        const auto& solved = accepted.front();
        const double pressure_scale = std::max({
            1.0,
            std::abs(solved.reference.block.pressure_pa),
            std::abs(solved.partner.block.pressure_pa),
        });
        result.diagnostics.pressure_residual_inf = std::abs(
            solved.reference.block.pressure_pa
            - solved.partner.block.pressure_pa
        ) / pressure_scale;
        for (std::size_t row = 0;
             row + 1 < solved.residual.size(); ++row) {
            result.diagnostics.chemical_potential_residual_inf =
                std::max(
                    result.diagnostics.chemical_potential_residual_inf,
                    std::abs(solved.residual[row])
                );
        }
        result.pressure_pa = 0.5 * (
            solved.reference.block.pressure_pa
            + solved.partner.block.pressure_pa
        );
        const PhaseSupport partner_support =
            problem.phases[partner_phase].support;
        result.phases.resize(2);
        result.phases[problem.reference_phase_index] = make_state(
            reference_eos,
            problem.phases[problem.reference_phase_index].role,
            PhaseSupport::AllComponents,
            solved.reference
        );
        result.phases[partner_phase] = make_state(
            partner_eos,
            problem.phases[partner_phase].role,
            partner_support,
            solved.partner
        );
        const auto& reference_state =
            result.phases[problem.reference_phase_index];
        const auto& partner_state = result.phases[partner_phase];
        double charge_scale = 0.0;
        for (std::size_t component = 0;
             component < reference_charges.size(); ++component) {
            charge_scale += std::abs(reference_charges[component])
                * reference_state.mole_fractions[component];
        }
        double partner_charge = 0.0;
        double partner_charge_scale = 0.0;
        for (std::size_t component = 0;
             component < partner_charges.size(); ++component) {
            partner_charge += partner_charges[component]
                * partner_state.mole_fractions[component];
            partner_charge_scale += std::abs(partner_charges[component])
                * partner_state.mole_fractions[component];
        }
        double reference_charge = 0.0;
        for (std::size_t component = 0;
             component < reference_charges.size(); ++component) {
            reference_charge += reference_charges[component]
                * reference_state.mole_fractions[component];
        }
        result.diagnostics.electroneutrality_inf = std::max(
            std::abs(reference_charge), std::abs(partner_charge)
        );
        const bool charge_certified =
            audit_held2_tolerance(
                kHeld2ChargeBalance, reference_charge, charge_scale
            ).passed
            && audit_held2_tolerance(
                kHeld2ChargeBalance,
                partner_charge,
                partner_charge_scale
            ).passed;
        result.diagnostics.phase_distinctness = std::abs(std::log(
            reference_state.molar_volume_m3_per_mol
            / partner_state.molar_volume_m3_per_mol
        ));
        const bool stable =
            reference_state.mechanical_class == "strict_stable"
            && partner_state.mechanical_class == "strict_stable";
        const bool ordered =
            reference_state.role == partner_state.role
            || (reference_state.role == PhaseRole::Liquid
                ? reference_state.molar_volume_m3_per_mol
                    < partner_state.molar_volume_m3_per_mol
                : reference_state.molar_volume_m3_per_mol
                    > partner_state.molar_volume_m3_per_mol);
        if (!stable || !ordered
            || result.diagnostics.phase_distinctness
                <= kDistinctTolerance) {
            result.diagnostics.eos_evaluation_count =
                static_cast<int>(eos_evaluation_count);
            result.phases.clear();
            result.diagnostics.solver_status = "converged";
            result.diagnostics.numerical_status = "certified";
            result.diagnostics.physical_status = "rejected";
            result.diagnostics.failure_stage = "phase_role_certificate";
            result.diagnostics.failure_reason =
                "coexistence_roles_not_thermodynamically_certified";
            return result;
        }
        if (!(result.pressure_pa > 0.0)
            || result.diagnostics.pressure_residual_inf > kPressureTolerance
            || result.diagnostics.chemical_potential_residual_inf
                > kPotentialTolerance
            || !charge_certified) {
            result.diagnostics.eos_evaluation_count =
                static_cast<int>(eos_evaluation_count);
            result.phases.clear();
            result.diagnostics.solver_status = "converged";
            result.diagnostics.numerical_status = "uncertified";
            result.diagnostics.physical_status = "rejected";
            result.diagnostics.failure_stage = "physical_certificate";
            result.diagnostics.failure_reason =
                "coexistence_pressure_potential_or_charge_not_certified";
            return result;
        }
        result.diagnostics.solver_status = "converged";
        result.diagnostics.numerical_status = "certified";
        result.diagnostics.physical_status = "accepted";
        if (!reference_eos.active_coordinates().empty()) {
            result.sensitivity = coexistence_sensitivity(
                reference_eos,
                partner_eos,
                problem,
                partner_to_reference,
                reference_coordinates,
                partner_coordinates,
                solved
            );
            eos_evaluation_count +=
                result.sensitivity->eos_evaluation_count;
        }
        if (!problem.reference_composition_directions.empty()) {
            result.reference_direction_sensitivity =
                coexistence_reference_direction_sensitivity(
                    reference_eos,
                    partner_eos,
                    problem,
                    partner_to_reference,
                    reference_coordinates,
                    partner_coordinates,
                    solved
                );
        }
        result.temperature_sensitivity = coexistence_temperature_sensitivity(
            reference_eos,
            partner_eos,
            problem,
            partner_to_reference,
            reference_coordinates,
            partner_coordinates,
            solved
        );
        eos_evaluation_count +=
            result.temperature_sensitivity->eos_evaluation_count;
        result.diagnostics.eos_evaluation_count =
            static_cast<int>(eos_evaluation_count);
    } catch (const std::exception& error) {
        result.diagnostics.eos_evaluation_count =
            static_cast<int>(eos_evaluation_count);
        result.phases.clear();
        result.diagnostics.solver_status = "failed";
        result.diagnostics.numerical_status = "not_certified";
        result.diagnostics.physical_status = "rejected";
        result.diagnostics.failure_stage = "local_equilibrium";
        result.diagnostics.failure_reason = error.what();
    }
    return result;
}

}  // namespace epcsaft_engine::equilibrium
