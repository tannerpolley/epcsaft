#include "eos.hpp"

#include "density.hpp"
#include "electrolyte_reference.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace epcsaft {
namespace {

struct NeutralReferenceBasis final {
    struct IonPair final {
        std::size_t cation;
        std::size_t anion;
        std::array<int, 2> stoichiometry;
    };

    std::vector<std::size_t> neutral_components;
    std::vector<IonPair> ion_pairs;
    std::vector<double> reference_composition;
    std::vector<double> rows;

    // EqExclude: input and domain validation
    [[nodiscard]] std::size_t row_count() const {
        return reference_composition.empty()
            ? 0u
            : rows.size() / reference_composition.size();
    }
};

// EqExclude: input and domain validation
bool full_row_rank(
    const std::vector<double>& values,
    std::size_t rows,
    std::size_t columns
) {
    if (rows == 0 || columns == 0 || values.size() != rows * columns) {
        return false;
    }
    std::vector<double> matrix(values);
    std::size_t rank = 0;
    for (std::size_t column = 0; column < columns && rank < rows; ++column) {
        std::size_t pivot = rank;
        while (pivot < rows
               && std::abs(matrix[pivot * columns + column]) <= 1.0e-12) {
            ++pivot;
        }
        if (pivot == rows) {
            continue;
        }
        for (std::size_t item = column; item < columns; ++item) {
            std::swap(
                matrix[rank * columns + item],
                matrix[pivot * columns + item]
            );
        }
        const double divisor = matrix[rank * columns + column];
        for (std::size_t item = column; item < columns; ++item) {
            matrix[rank * columns + item] /= divisor;
        }
        for (std::size_t row = 0; row < rows; ++row) {
            if (row == rank) {
                continue;
            }
            const double factor = matrix[row * columns + column];
            for (std::size_t item = column; item < columns; ++item) {
                matrix[row * columns + item] -=
                    factor * matrix[rank * columns + item];
            }
        }
        ++rank;
    }
    return rank == rows;
}

// EqID: neutral_reference_stoichiometric_basis
NeutralReferenceBasis make_neutral_reference_basis(
    const std::vector<std::string>& component_ids,
    const std::vector<int>& charges,
    const std::vector<double>& reference_composition,
    const std::vector<NeutralReferenceBasis::IonPair>& ion_pairs
) {
    NeutralReferenceBasis result;
    const std::size_t components = component_ids.size();
    if (components < 3 || components != charges.size()
        || reference_composition.size() != components
        || std::any_of(
            component_ids.begin(), component_ids.end(), [](const auto& value) {
                return value.empty();
            }
        )) {
        return result;
    }
    std::vector<std::size_t> order(components);
    std::iota(order.begin(), order.end(), 0u);
    std::sort(order.begin(), order.end(), [&](std::size_t left, std::size_t right) {
        return component_ids[left] < component_ids[right];
    });
    for (std::size_t index = 1; index < order.size(); ++index) {
        if (component_ids[order[index - 1]] == component_ids[order[index]]) {
            return {};
        }
    }
    for (std::size_t component : order) {
        if (charges[component] == 0) {
            result.neutral_components.push_back(component);
        }
    }
    if (result.neutral_components.empty() || ion_pairs.empty()) {
        return {};
    }
    result.reference_composition = reference_composition;
    auto append_row = [&](std::size_t left, double left_value,
                          std::size_t right = 0u, double right_value = 0.0) {
        result.rows.resize(result.rows.size() + components, 0.0);
        const std::size_t offset = result.rows.size() - components;
        result.rows[offset + left] = left_value;
        if (right_value != 0.0) {
            result.rows[offset + right] = right_value;
        }
    };
    for (std::size_t component : result.neutral_components) {
        append_row(component, 1.0);
    }
    for (const auto& pair : ion_pairs) {
        if (pair.cation >= components || pair.anion >= components
            || charges[pair.cation] <= 0 || charges[pair.anion] >= 0
            || pair.stoichiometry[0] <= 0 || pair.stoichiometry[1] <= 0) {
            return {};
        }
        const std::int64_t charge_balance =
            static_cast<std::int64_t>(pair.stoichiometry[0])
                * static_cast<std::int64_t>(charges[pair.cation])
            + static_cast<std::int64_t>(pair.stoichiometry[1])
                * static_cast<std::int64_t>(charges[pair.anion]);
        if (charge_balance != 0) {
            return {};
        }
        result.ion_pairs.push_back(pair);
        append_row(
            pair.cation,
            static_cast<double>(pair.stoichiometry[0]),
            pair.anion,
            static_cast<double>(pair.stoichiometry[1])
        );
    }
    const std::size_t expected_rows = components - 1u;
    if (result.row_count() != expected_rows
        || !full_row_rank(result.rows, expected_rows, components)) {
        return {};
    }
    for (std::size_t row = 0; row < expected_rows; ++row) {
        double net_charge = 0.0;
        for (std::size_t component = 0; component < components; ++component) {
            net_charge += result.rows[row * components + component]
                * static_cast<double>(charges[component]);
        }
        if (net_charge != 0.0) {
            return {};
        }
    }
    return result;
}

// EqExclude: input and domain validation
bool sha256_fingerprint_valid(const std::string& fingerprint) {
    return fingerprint.size() == 71
        && fingerprint.rfind("sha256:", 0) == 0
        && std::all_of(
            fingerprint.begin() + 7,
            fingerprint.end(),
            [](unsigned char character) {
                return std::isxdigit(character) != 0;
            }
        );
}

}  // namespace

// EqID: declared_solvent_source_reference_transfer
NativeSourceReferenceResult NativeModel::evaluate_source_reference_transfer(
    double temperature_k,
    double pressure_pa,
    const NativeSourceReferenceRequest& request
) const {
    NativeSourceReferenceResult result;
    const auto fail = [&](int status, const std::string& message) {
        result.status = status;
        result.error = message;
        return result;
    };
    if (!supports_electrolyte_transport()
        || eos_metadata().topology_fingerprint.empty()
        || !std::isfinite(eos_metadata().temperature_min_k)
        || !std::isfinite(eos_metadata().temperature_max_k)
        || !std::isfinite(eos_metadata().pressure_min_pa)
        || !std::isfinite(eos_metadata().pressure_max_pa)) {
        return fail(
            SOURCE_REFERENCE_MISSING_MODEL,
            "source-reference model identity is incompatible"
        );
    }
    const bool candidate = !eos_metadata().candidate_qualification.empty();
    if (candidate && (
            eos_metadata().candidate_qualification != "candidate_extrapolation"
            || !sha256_fingerprint_valid(eos_metadata().candidate_source_fingerprint)
            || !std::isfinite(eos_metadata().source_temperature_min_k)
            || !std::isfinite(eos_metadata().source_temperature_max_k)
            || !std::isfinite(eos_metadata().source_pressure_min_pa)
            || !std::isfinite(eos_metadata().source_pressure_max_pa))) {
        return fail(
            SOURCE_REFERENCE_MISSING_MODEL,
            "candidate source-reference identity is incomplete"
        );
    }
    if (!std::isfinite(temperature_k) || temperature_k <= 0.0
        || !std::isfinite(pressure_pa) || pressure_pa <= 0.0) {
        return fail(
            SOURCE_REFERENCE_NONFINITE,
            "source-reference temperature and pressure must be positive and finite"
        );
    }
    if (temperature_k < eos_metadata().temperature_min_k
        || temperature_k > eos_metadata().temperature_max_k
        || pressure_pa < eos_metadata().pressure_min_pa
        || pressure_pa > eos_metadata().pressure_max_pa) {
        return fail(
            SOURCE_REFERENCE_DOMAIN_REJECTION,
            "source-reference trial state is outside the EOS domain"
        );
    }
    const std::size_t component_count = this->component_count();
    const std::size_t row_count = component_count - 1u;
    if (request.component_ids.size() != component_count
        || request.source_solvent_mole_fractions.size() != component_count
        || request.ion_pair_indices.empty()
        || request.ion_pair_indices.size()
            != request.ion_pair_stoichiometry.size()) {
        return fail(
            SOURCE_REFERENCE_COMPONENT_MISMATCH,
            "source-reference component dimensions do not match"
        );
    }
    if (request.component_ids != component_ids()) {
        return fail(
            SOURCE_REFERENCE_COMPONENT_MISMATCH,
            "source-reference component order does not match"
        );
    }
    const bool molality_reference =
        request.activity_convention_id
            == "molality-infinite-dilution-v1"
        && request.reference_convention == "molality-infinite-dilution";
    const bool mixed_mole_fraction_reference =
        request.activity_convention_id
            == "mole-fraction-mixed-reference-v1"
        && request.reference_convention == "mole-fraction-mixed-reference";
    if (request.reference_state_id.empty()
        || (!molality_reference && !mixed_mole_fraction_reference)
        || request.phase != "liquid"
        || !sha256_fingerprint_valid(request.reference_state_fingerprint)
        || !std::isfinite(request.standard_molality_mol_per_kg)
        || request.standard_molality_mol_per_kg <= 0.0
        || !std::isfinite(request.source_reference_pressure_pa)
        || request.source_reference_pressure_pa <= 0.0) {
        return fail(
            SOURCE_REFERENCE_UNSUPPORTED,
            "source-reference declaration or activity convention is unsupported"
        );
    }
    constexpr std::uint32_t pressure_derivative = 2u;
    constexpr std::uint32_t parameter_derivative = 8u;
    constexpr std::uint32_t supported_derivatives =
        pressure_derivative | parameter_derivative;
    const bool pressure_derivatives_requested =
        (request.required_derivatives & pressure_derivative) != 0u;
    const bool parameter_derivatives_requested =
        (request.required_derivatives & parameter_derivative) != 0u;
    if ((request.required_derivatives & ~supported_derivatives) != 0u
        || (parameter_derivatives_requested
            && active_coordinates().empty())) {
        return fail(
            SOURCE_REFERENCE_DERIVATIVE_UNAVAILABLE,
            "required source-reference derivative is unavailable"
        );
    }
    const std::size_t active_parameter_count =
        active_coordinates().size();
    try {
        std::vector<double> source_composition(
            request.source_solvent_mole_fractions
        );
        double sum = 0.0;
        for (std::size_t component = 0; component < component_count; ++component) {
            const double fraction = source_composition[component];
            if (!std::isfinite(fraction) || fraction < 0.0
                || (charge_numbers()[component] != 0
                    && fraction != 0.0)) {
                return fail(
                    SOURCE_REFERENCE_INVALID_SOURCE,
                    "source composition must be finite, salt-free, and nonnegative"
                );
            }
            sum += fraction;
        }
        if (std::abs(sum - 1.0) > 1.0e-12) {
            return fail(
                SOURCE_REFERENCE_INVALID_SOURCE,
                "source composition must sum to one"
            );
        }
        if (request.pure_component_indices.size()
            != request.pure_component_phases.size()) {
            return fail(
                SOURCE_REFERENCE_INVALID_SOURCE,
                "pure-component reference identities or phases are unavailable"
            );
        }
        std::vector<bool> pure_component(component_count, false);
        std::vector<std::string> pure_component_phase(component_count);
        for (std::size_t index = 0;
             index < request.pure_component_indices.size();
             ++index) {
            const std::size_t component = request.pure_component_indices[index];
            const std::string& phase = request.pure_component_phases[index];
            if (component >= component_count
                || charge_numbers()[component] != 0
                || pure_component[component]
                || (phase != "liquid" && phase != "vapor")) {
                return fail(
                    SOURCE_REFERENCE_INVALID_SOURCE,
                    "pure-component references must be unique neutral components "
                    "with liquid or vapor phases"
                );
            }
            pure_component[component] = true;
            pure_component_phase[component] = phase;
        }
        if (molality_reference && !request.pure_component_indices.empty()) {
            return fail(
                SOURCE_REFERENCE_INVALID_SOURCE,
                "molality reference does not admit pure-component overrides"
            );
        }
        std::vector<NeutralReferenceBasis::IonPair> ion_pairs;
        ion_pairs.reserve(request.ion_pair_indices.size());
        for (std::size_t index = 0;
             index < request.ion_pair_indices.size();
             ++index) {
            ion_pairs.push_back({
                request.ion_pair_indices[index][0],
                request.ion_pair_indices[index][1],
                request.ion_pair_stoichiometry[index],
            });
        }
        const auto reference_basis = make_neutral_reference_basis(
            component_ids(),
            charge_numbers(),
            source_composition,
            ion_pairs
        );
        if (reference_basis.row_count() != row_count) {
            return fail(
                SOURCE_REFERENCE_INVALID_SOURCE,
                "source-reference ion pairs do not define a full neutral basis"
            );
        }
        const auto& molar_masses = molar_masses_kg_per_mol();
        const double solvent_molar_mass = std::inner_product(
            source_composition.begin(),
            source_composition.end(),
            molar_masses.begin(),
            0.0
        );
        if (!std::isfinite(solvent_molar_mass)
            || solvent_molar_mass <= 0.0) {
            return fail(
                SOURCE_REFERENCE_INVALID_SOURCE,
                "source solvent molar mass is unavailable"
            );
        }
        std::vector<double> fugacity(row_count, 0.0);
        std::vector<double> activity(row_count, 0.0);
        std::vector<double> transfer(row_count, 0.0);
        std::vector<double> fugacity_pressure_derivatives(
            pressure_derivatives_requested ? row_count : 0u, 0.0
        );
        std::vector<double> fugacity_parameter_derivatives(
            parameter_derivatives_requested
                ? row_count * active_parameter_count
                : 0u,
            0.0
        );
        const std::size_t derivative_stride = 1u
            + (pressure_derivatives_requested ? 1u : 0u)
            + (parameter_derivatives_requested ? active_parameter_count : 0u);
        std::vector<std::size_t> neutral_solutes;
        for (std::size_t component : reference_basis.neutral_components) {
            if (source_composition[component] == 0.0
                && !pure_component[component]) {
                neutral_solutes.push_back(component);
            }
        }
        std::vector<double> neutral_trace_stoichiometry(component_count, 0.0);
        for (std::size_t component : neutral_solutes) {
            neutral_trace_stoichiometry[component] = 1.0;
        }
        const auto reference_composition = [&] (
            const std::vector<double>& trace_stoichiometry,
            double molality
        ) {
            const double trace_formula_moles = std::accumulate(
                trace_stoichiometry.begin(),
                trace_stoichiometry.end(),
                0.0
            );
            const double solvent_moles = 1.0 / solvent_molar_mass;
            const double total_moles = solvent_moles
                + molality * trace_formula_moles;
            std::vector<double> composition(component_count, 0.0);
            for (std::size_t component = 0;
                 component < component_count;
                 ++component) {
                composition[component] =
                    source_composition[component] * solvent_moles / total_moles
                    + molality * trace_stoichiometry[component] / total_moles;
            }
            return composition;
        };
        std::optional<epcsaft::NativeResult> density_anchor;
        if (pressure_derivatives_requested
            || parameter_derivatives_requested) {
            const std::vector<double> density_anchor_composition =
                reference_composition(
                    neutral_trace_stoichiometry,
                    detail::REFERENCE_MOLALITIES[2]
                );
            const NativeResult anchor = evaluate_pressure(
                temperature_k,
                pressure_pa,
                density_anchor_composition,
                "liquid"
            );
            result.evaluations.state += anchor.eos_evaluation_count;
            density_anchor = anchor;
        }
        const auto phase_fugacity_derivatives =
            [&](const std::vector<double>& composition,
                const char* phase_branch,
                bool use_density_anchor) {
                const auto state = density_anchor.has_value()
                        && use_density_anchor
                    ? evaluate_pressure_continuation(
                        temperature_k,
                        pressure_pa,
                        composition,
                        density_anchor->molar_density,
                        phase_branch,
                        density_anchor->density_diagnostics
                            ->selected_stable_root_index,
                        density_anchor->density_diagnostics
                            ->certified_origin_molar_density
                    )
                    : evaluate_pressure(
                        temperature_k,
                        pressure_pa,
                        composition,
                        phase_branch
                    );
                result.evaluations.state += state.eos_evaluation_count;
                if (!state.log_fugacity_coefficient.has_value()) {
                    throw std::runtime_error(
                        "source-reference fugacity is unavailable"
                    );
                }
                std::vector<double> values(
                    component_count * derivative_stride, 0.0
                );
                for (std::size_t component = 0;
                     component < component_count;
                     ++component) {
                    values[component * derivative_stride] =
                        (*state.log_fugacity_coefficient)[component];
                }
                if (!pressure_derivatives_requested
                    && !parameter_derivatives_requested) {
                    return values;
                }
                const auto absent_coordinate =
                    std::numeric_limits<std::size_t>::max();
                std::vector<std::size_t> support_coordinates(
                    component_count, absent_coordinate
                );
                std::vector<std::size_t> supported_component_indices;
                supported_component_indices.reserve(component_count);
                std::size_t support_count = 0u;
                for (std::size_t component = 0;
                     component < component_count;
                     ++component) {
                    if (composition[component] > 0.0) {
                        support_coordinates[component] = support_count++;
                        supported_component_indices.push_back(component);
                    }
                    std::fill(
                        values.begin()
                            + static_cast<std::ptrdiff_t>(
                                component * derivative_stride + 1u
                            ),
                        values.begin()
                            + static_cast<std::ptrdiff_t>(
                                (component + 1u) * derivative_stride
                            ),
                        std::numeric_limits<double>::quiet_NaN()
                    );
                }
                const double volume_m3 = 1.0 / state.molar_density;
                const NativeMixturePhaseBlock phase =
                    evaluate_supported_phase_derivatives(
                        temperature_k,
                        composition,
                        volume_m3,
                        eos_metadata().ion_mole_fraction_max,
                        active_coordinates(),
                        supported_component_indices
                    );
                ++result.evaluations.state;
                const std::size_t state_count = support_count + 1u;
                const std::size_t width =
                    state_count + active_parameter_count;
                if (phase.hessian.size() != width * width) {
                    throw std::runtime_error(
                        "source-reference derivative tensor dimensions changed"
                    );
                }
                const std::size_t volume_coordinate = support_count;
                const double pressure_volume_derivative =
                    -detail::GAS_CONSTANT * temperature_k
                    * phase.hessian[
                        volume_coordinate * width + volume_coordinate
                    ];
                if (!std::isfinite(pressure_volume_derivative)
                    || pressure_volume_derivative == 0.0) {
                    throw std::runtime_error(
                        "source-reference pressure closure is singular"
                    );
                }
                for (std::size_t component = 0;
                     component < component_count;
                     ++component) {
                    const std::size_t component_coordinate =
                        support_coordinates[component];
                    if (component_coordinate == absent_coordinate) {
                        continue;
                    }
                    std::size_t offset = component * derivative_stride + 1u;
                    const double chemical_volume = phase.hessian[
                        component_coordinate * width + volume_coordinate
                    ];
                    if (pressure_derivatives_requested) {
                        values[offset++] = chemical_volume
                            / pressure_volume_derivative
                            - 1.0 / pressure_pa;
                    }
                    if (parameter_derivatives_requested) {
                        for (std::size_t parameter = 0;
                             parameter < active_parameter_count;
                             ++parameter) {
                            const std::size_t parameter_coordinate =
                                state_count + parameter;
                            values[offset + parameter] = phase.hessian[
                                component_coordinate * width
                                    + parameter_coordinate
                            ] - chemical_volume * phase.hessian[
                                volume_coordinate * width
                                    + parameter_coordinate
                            ] / phase.hessian[
                                volume_coordinate * width + volume_coordinate
                            ];
                        }
                    }
                }
                return values;
            };
        const auto selected_derivatives = [&] (
            const std::vector<double>& values,
            const std::vector<std::size_t>& components
        ) {
            std::vector<double> selected;
            selected.reserve(components.size() * derivative_stride);
            for (std::size_t component : components) {
                const auto first = values.begin()
                    + static_cast<std::ptrdiff_t>(
                        component * derivative_stride
                    );
                if (!std::all_of(
                        first,
                        first + static_cast<std::ptrdiff_t>(derivative_stride),
                        [](double value) { return std::isfinite(value); }
                    )) {
                    throw std::runtime_error(
                        "source-reference derivative component is outside "
                        "the positive limiting support"
                    );
                }
                selected.insert(
                    selected.end(),
                    first,
                    first + static_cast<std::ptrdiff_t>(derivative_stride)
                );
            }
            return selected;
        };
        const auto assign_derivatives = [&] (
            std::size_t target_row,
            const std::vector<double>& values,
            std::size_t selected_component
        ) {
            std::size_t offset = selected_component * derivative_stride;
            fugacity[target_row] = values[offset++];
            if (pressure_derivatives_requested) {
                fugacity_pressure_derivatives[target_row] = values[offset++];
            }
            if (parameter_derivatives_requested) {
                for (std::size_t parameter = 0;
                     parameter < active_parameter_count;
                     ++parameter) {
                    fugacity_parameter_derivatives[
                        target_row * active_parameter_count + parameter
                    ] = values[offset + parameter];
                }
            }
        };
        double reference_molality =
            detail::REFERENCE_MOLALITIES.back();
        double convergence_error = 0.0;
        std::vector<double> neutral_values(
            reference_basis.neutral_components.size() * derivative_stride,
            0.0
        );
        std::vector<std::size_t> non_pure_neutral_components;
        for (std::size_t component : reference_basis.neutral_components) {
            if (!pure_component[component]) {
                non_pure_neutral_components.push_back(component);
            }
        }
        const auto copy_neutral_values = [&] (
            const std::vector<std::size_t>& selected_components,
            const std::vector<double>& values
        ) {
            for (std::size_t selected = 0;
                 selected < selected_components.size();
                 ++selected) {
                const auto position = std::find(
                    reference_basis.neutral_components.begin(),
                    reference_basis.neutral_components.end(),
                    selected_components[selected]
                );
                if (position == reference_basis.neutral_components.end()) {
                    throw std::runtime_error(
                        "source-reference neutral basis changed"
                    );
                }
                const std::size_t target = static_cast<std::size_t>(
                    std::distance(
                        reference_basis.neutral_components.begin(),
                        position
                    )
                );
                std::copy_n(
                    values.begin()
                        + static_cast<std::ptrdiff_t>(
                            selected * derivative_stride
                        ),
                    derivative_stride,
                    neutral_values.begin()
                        + static_cast<std::ptrdiff_t>(
                            target * derivative_stride
                        )
                );
            }
        };
        const auto reference_limit_derivatives = [&] (
            const std::vector<std::size_t>& selected_components,
            const std::vector<double>& trace_stoichiometry
        ) {
            const auto evaluate_limit = [&](double molality) {
                const std::vector<double> composition =
                    reference_composition(trace_stoichiometry, molality);
                return selected_derivatives(
                    phase_fugacity_derivatives(composition, "liquid", true),
                    selected_components
                );
            };
            return pressure_derivatives_requested
                    || parameter_derivatives_requested
                ? detail::reference_limit_terminal_pair(evaluate_limit)
                : detail::reference_limit_sequence(evaluate_limit);
        };
        const double neutral_trace_formula_moles = std::accumulate(
            neutral_trace_stoichiometry.begin(),
            neutral_trace_stoichiometry.end(),
            0.0
        );
        if (neutral_trace_formula_moles == 0.0) {
            copy_neutral_values(
                non_pure_neutral_components,
                selected_derivatives(
                    phase_fugacity_derivatives(
                        source_composition, "liquid", true
                    ),
                    non_pure_neutral_components
                )
            );
        } else {
            const auto sequence = reference_limit_derivatives(
                non_pure_neutral_components,
                neutral_trace_stoichiometry
            );
            copy_neutral_values(non_pure_neutral_components, sequence.values);
            reference_molality = sequence.molality;
            convergence_error = sequence.convergence_error;
        }
        for (std::size_t neutral = 0;
             neutral < reference_basis.neutral_components.size();
             ++neutral) {
            const std::size_t component =
                reference_basis.neutral_components[neutral];
            if (!pure_component[component]) {
                continue;
            }
            std::vector<double> composition(component_count, 0.0);
            composition[component] = 1.0;
            const std::vector<double> values = selected_derivatives(
                phase_fugacity_derivatives(
                    composition,
                    pure_component_phase[component].c_str(),
                    false
                ),
                {component}
            );
            std::copy(
                values.begin(),
                values.end(),
                neutral_values.begin()
                    + static_cast<std::ptrdiff_t>(neutral * derivative_stride)
            );
        }
        const double solute_scale = std::log(
            solvent_molar_mass * request.standard_molality_mol_per_kg
        );
        std::size_t row = 0;
        for (std::size_t component : reference_basis.neutral_components) {
            assign_derivatives(row, neutral_values, row);
            if (mixed_mole_fraction_reference) {
                activity[row] = 0.0;
            } else if (source_composition[component] > 0.0) {
                activity[row] = std::log(source_composition[component]);
            } else {
                activity[row] = solute_scale;
            }
            transfer[row] = fugacity[row] + activity[row];
            ++row;
        }
        const auto append_pair = [&](const NeutralReferenceBasis::IonPair& pair) {
            const std::size_t cation = pair.cation;
            const std::size_t anion = pair.anion;
            const auto& stoichiometry = pair.stoichiometry;
            auto pair_trace_stoichiometry = neutral_trace_stoichiometry;
            pair_trace_stoichiometry[cation] +=
                static_cast<double>(stoichiometry[0]);
            pair_trace_stoichiometry[anion] +=
                static_cast<double>(stoichiometry[1]);
            const auto sequence = reference_limit_derivatives(
                {cation, anion}, pair_trace_stoichiometry
            );
            for (std::size_t derivative = 0;
                 derivative < derivative_stride;
                 ++derivative) {
                const double contraction =
                    stoichiometry[0] * sequence.values[derivative]
                    + stoichiometry[1]
                        * sequence.values[derivative_stride + derivative];
                if (derivative == 0u) {
                    fugacity[row] = contraction;
                } else if (pressure_derivatives_requested
                           && derivative == 1u) {
                    fugacity_pressure_derivatives[row] = contraction;
                } else {
                    const std::size_t parameter = derivative - 1u
                        - (pressure_derivatives_requested ? 1u : 0u);
                    fugacity_parameter_derivatives[
                        row * active_parameter_count + parameter
                    ] = contraction;
                }
            }
            activity[row] = mixed_mole_fraction_reference
                ? 0.0
                : (static_cast<double>(stoichiometry[0])
                   + static_cast<double>(stoichiometry[1])) * solute_scale;
            transfer[row] = fugacity[row] + activity[row];
            reference_molality = sequence.molality;
            convergence_error = std::max(
                convergence_error, sequence.convergence_error
            );
            ++row;
        };
        for (const auto& pair : reference_basis.ion_pairs) {
            append_pair(pair);
        }
        const auto finite = [](const std::vector<double>& values) {
            return std::all_of(values.begin(), values.end(), [](double value) {
                return std::isfinite(value);
            });
        };
        if (row != row_count || !finite(fugacity) || !finite(activity)
            || !finite(transfer) || !finite(fugacity_pressure_derivatives)
            || !finite(fugacity_parameter_derivatives)
            || !std::isfinite(convergence_error)) {
            return fail(
                SOURCE_REFERENCE_NONFINITE,
                "source-reference transfer is nonfinite"
            );
        }
        result.component_count = component_count;
        result.neutral_basis_row_count = row_count;
        result.neutral_basis = reference_basis.rows;
        result.log_fugacity_contractions = std::move(fugacity);
        result.activity_scale_log_contractions = std::move(activity);
        result.transfer_log_contractions = std::move(transfer);
        result.source_composition = std::move(source_composition);
        result.derivative_availability = request.required_derivatives;
        result.pressure_derivatives_per_pa =
            std::move(fugacity_pressure_derivatives);
        result.parameter_derivatives =
            std::move(fugacity_parameter_derivatives);
        result.temperature_k = temperature_k;
        result.pressure_pa = pressure_pa;
        result.source_reference_pressure_pa =
            request.source_reference_pressure_pa;
        result.standard_molality_mol_per_kg =
            request.standard_molality_mol_per_kg;
        result.solvent_molar_mass_kg_per_mol = solvent_molar_mass;
        result.reference_limit_molality_mol_per_kg = reference_molality;
        result.reference_convergence_error = convergence_error;
        result.source_temperature_min_k = eos_metadata().source_temperature_min_k;
        result.source_temperature_max_k = eos_metadata().source_temperature_max_k;
        result.source_pressure_min_pa = eos_metadata().source_pressure_min_pa;
        result.source_pressure_max_pa = eos_metadata().source_pressure_max_pa;
        result.candidate_temperature_min_k = eos_metadata().temperature_min_k;
        result.candidate_temperature_max_k = eos_metadata().temperature_max_k;
        result.candidate_pressure_min_pa = eos_metadata().pressure_min_pa;
        result.candidate_pressure_max_pa = eos_metadata().pressure_max_pa;
        result.candidate_qualification = eos_metadata().candidate_qualification;
        result.candidate_source_fingerprint =
            eos_metadata().candidate_source_fingerprint;
        result.parameter_fingerprint = eos_metadata().parameter_fingerprint;
        result.topology_fingerprint = eos_metadata().topology_fingerprint;
        result.component_order_fingerprint =
            eos_metadata().component_order_fingerprint;
        result.reference_state_fingerprint =
            request.reference_state_fingerprint;
        result.domain_fingerprint = eos_metadata().domain_fingerprint;
        result.artifact_fingerprint = eos_metadata().artifact_fingerprint;
        result.helmholtz_basis_id = EOS_HELMHOLTZ_BASIS_ID;
        result.status = SOURCE_REFERENCE_OK;
        return result;
    } catch (const detail::DensityDomainError& error) {
        return fail(SOURCE_REFERENCE_DOMAIN_REJECTION, error.what());
    } catch (const std::invalid_argument& error) {
        return fail(SOURCE_REFERENCE_INVALID_SOURCE, error.what());
    } catch (const std::runtime_error& error) {
        const std::string message(error.what());
        const bool nonfinite = message.find("nonfinite") != std::string::npos;
        return fail(
            nonfinite
                ? SOURCE_REFERENCE_NONFINITE
                : (message.find("limit") == std::string::npos
                    ? SOURCE_REFERENCE_EVALUATION_ERROR
                    : SOURCE_REFERENCE_LIMIT_FAILURE),
            error.what()
        );
    } catch (const std::exception& error) {
        return fail(SOURCE_REFERENCE_EVALUATION_ERROR, error.what());
    } catch (...) {
        return fail(
            SOURCE_REFERENCE_EVALUATION_ERROR,
            "unknown source-reference transfer failure"
        );
    }
}

}  // namespace epcsaft
